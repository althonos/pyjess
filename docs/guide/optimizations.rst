Optimizations
=============

.. currentmodule:: pyjess

While PyJess started as a Cython wrapper of Jess, it also 
contains several optimizations to make the code write better
while maintaining consistency with the original Jess code.
Some of these optimizations are described below, as well as
the version where they were introduced.


Residue name index
------------------

.. versionadded:: 0.5.1

Given a `Molecule`, matching a new `Template` to that `Molecule`
requires generating a set of candidate atoms for each `TemplateAtom`
of the template according to its `~TemplateAtom.match_mode`. 

In the original Jess code, this is done by a full-scan on each 
`Atom` of the `Molecule` for each `TemplateAtom` at index ``k``:

.. code:: c

    Atom* A;
    Molecule* M;
    Template* T;
    CandidateSet* S;
    int n, i;

    n=Molecule_count(M);
    for(i=0; m<n; m++)
    {
        A=(Atom*) Molecule_atom(M,i);
        if(T->match(T,k,A))
        {
            S->atom[S->count]=A;
            S->count++;
        }
    }

However, the most common match modes require the candidate `Atom` to 
have a `~Atom.residue_name` equal to one of the `TemplateAtom` 
`~TemplateAtom.residue_names`. To exploit this, and avoid a full-scan,
we can create an index grouping the atoms of a `Molecule` by their
`~Atom.residue_name`, and only iterate on atoms with the right 
`~Atom.residue_name` when building the candidates for a `TemplateAtom`:

.. code:: c

    Atom* A;
    Molecule* M;
    Template* T;
    CandidateSet* S;
    Atom** atoms;

    for(name=T->residueNames(T,k);name!=NULL;name++)
    {
        for(atoms=Molecule_atoms(M,resname);*atoms!=NULL;atoms++)
        {
            A=*atoms;
            if(T->match(T,k,A))
            {
                S->atom[S->count]=A;
                S->count++;
            }
        }
    }

By doing so we greatly reduce the number of calls to ``T->match``, 
at the cost of computing an index when a new `Molecule` is created.
This is grealy beneficial for a large number of templates, with an 
additional :math:`O(n)` memory requirement. 


k-d tree generation
-------------------

.. versionadded:: 0.6.0

Jess uses a `k-d tree <https://en.wikipedia.org/wiki/K-d_tree>`_ data 
structure to partition the molecule into geometric regions and speed-up the 
retrieval of atoms in each region. To generate a :math:`k`-d tree from 
a list of `Atom`, the atoms are recursively partitioned on a single dimension
around a pivot value, usually the median coordinate for that dimension.

The original Jess code uses ``qsort`` (the `QuickSort <https://en.wikipedia.org/wiki/Quicksort>`_
implementation of the C standard library) to first sort the atoms on a 
single dimension, and then takes the middle point:

.. code:: c

    KdTreeNode* N;
    int* idx;
    int n;

    qsort(idx,n,sizeof(int),KdTree_compare);

    split = n/2;
    N->index=idx[split-1];
    N->type=type;

    type = (type+1)%dim;
    N->left = KdTreeNode_create(idx,split,type,u,dim);
    N->right = KdTreeNode_create(&idx[split],n-split,type,u,dim); 


While implemented very efficiently, QuickSort has an average runtime
complexity of :math:`O(nlog(n))`. Given that the algorithm is only used
here to search for the median, and that the sort order is actually irrelevant, 
we replaced it with the `QuickSelect <https://en.wikipedia.org/wiki/Quickselect>`_ algorithm,
which can retrieve the median with an average runtime complexity of :math:`O(n)`.


Approximate annulus intersection
--------------------------------

.. versionadded:: 0.6.0

During the search for matches in a ``Scanner``, Jess takes into account the 
flexibility of the `Template` by modelling the distance constraints to each 
`TemplateAtom` using an `annulus <https://en.wikipedia.org/wiki/Annulus_(mathematics)>`_. 
Then, it queries the :math:`k`-d tree created on the candidate atoms to find
candidate atoms that are included in this annulus.

In the original Jess, the traversal of the :math:`k`-d tree is done by computing 
for each internal node of the tree whether the box formed by that node intersect
the query annulus, and for each leaf whether they are contained in that annulus.

Computing the intersection between a box and an annulus requires computing 
`Euclidean distances <https://en.wikipedia.org/wiki/Euclidean_distance>`_, 
and therefore products between real numbers, an operation that is among the 
slowest even on modern CPUs:

.. code:: c

    Annulus* A;
    double minBox[d];
    double maxBox[d];

    double minSum=0.0;
    double maxSum=0.0;
    for(i=0; i<d; i++)
    {
        double t1 = A->centre[i]-minBox[i];
        double t2 = A->centre[i]-maxBox[i];
        t1 *= t1;
        t2 *= t2;

        if(minBox[i]>A->centre[i] || maxBox[i]<A->centre[i])
            minSum += min(t1,t2);

        maxSum += max(t1,t2);
    }

    bool intersects = minSum>(A->outer*A->outer) || maxSum<(A->inner*A->inner);

To speed-up the querying, we approximate the query annulus as a bounding
box, and instead compute the intersection to the :math:`k`-d tree box 
to the bounding box, which only requires comparisons:

.. code:: c

    Annulus* A;
    double minBox[d];
    double maxBox[d];

    bool intersects = true;
    for(i=0; i<d; i++)
    {
        double dmin = A->centre[i]-A->outer;
        double dmax = A->centre[i]+A->outer;
        if( !( dmin <= maxBox[i] && minBox[i] <= dmax ) )
            intersects = false;
    }

As this is an approximation, it may wrongly return ``true`` on (literal) corner 
cases, i.e. when the intersection happens in a corner of the bounding box around
the annulus. However, given that the :math:`k`-d tree later checks that the 
points are actually included in the annulus for the leaf nodes, 
using this implementation will not generate false positives.


Reduced backtracking
--------------------

.. versionadded:: 0.6.0

The ``Scanner`` in Jess iterates on the `Template` atoms and then aligns the
`Molecule` atoms iteratively. When no more candidate for a `TemplateAtom` can
be found in a `Molecule`, it backtracks to the previous `TemplateAtom`, and 
continues this iteration. 

Since every backtracking event triggers the querying of the :math:`k`-d tree, 
we want to minimize backtracking as much as possible. To do so, we compute 
an iteration order over the `Template` atoms that minimizes the amount of 
backtracking required to explore all paths. This can be done quickly by 
sorting the `TemplateAtom` using the number of `Atom` in the corresponding 
``CandidateSet``. 

Empirically, this approach reduced the amount of :math:`k`-d tree queries by a 
factor of 10. It is nevertheless unclear whether an optimal path can be 
further identified and pre-computed from the candidate `Atom`, possibly by 
filtering distant atoms first.

.. warning:: 

    Out all the optimizations in this section, this one is the only one 
    to introduce *slight* behavioral changes in PyJess compared to Jess.
    Since the order in which candidate atoms are matched have changed, 
    the order in which a PyJess `Query` yields `Hit` objects differs
    to that in which Jess reports a match. This only affects the *order*
    though: all matches are still returned, and are 1-to-1 identical!
    This can be disabled by running with ``Jess.query(..., reorder=False)``
    to use the original matching order, at the cost of a longer runtime.



Type concretization
-------------------

.. versionadded:: 0.6.0

Jess uses generic code to support multiple types of spatial regions
using function pointers to implement `virtual method tables <https://en.wikipedia.org/wiki/Virtual_method_table>`_.
While elegant and functional, the code never really makes uses of 
the genericity, and therefore the algorithm can be specialized to 
the appropriate `Region` concrete type (either `Join` or `Annulus`)
to remove the overhead of calling the function pointers. 

In addition, most of the geometric code is inlined, as it is called in hot paths
where the compiler can apply auto-vectorization.


Dimension concretization
------------------------

.. versionadded:: 0.6.0

The `Annulus` type is made to be generic over the number of dimensions,
but in practice it is only used for 3-dimensions. We hardcoded the dimensions
to encourage the compiler to unroll loops over the `Annulus` dimensions where
applicable, and to use constant-size arrays rather than dynamic allocation when
applicable.


Scanner memory recycling
------------------------

.. versionadded:: 0.6.0

The original Jess code performs a lot of allocation/deallocations in hot paths,
as it creates a new `Scanner` which allocates memory for each `Template` / `Molecule`
pair to match. Effectively, most of these buffers can actually be reused 
across `Templates` for a given `Molecule`, provided sufficient bookkeeping. Our
implementation keeps allocation to a minimum across an entire `Query`.
