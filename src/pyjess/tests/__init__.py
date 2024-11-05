from . import (
    test_atom,
    test_jess,
    test_molecule,
    test_template_atom,
    test_template
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_atom))
    suite.addTests(loader.loadTestsFromModule(test_jess))
    suite.addTests(loader.loadTestsFromModule(test_molecule))
    suite.addTests(loader.loadTestsFromModule(test_template_atom))
    suite.addTests(loader.loadTestsFromModule(test_template))
    # test_doctest.load_tests(loader, suite, pattern)
    return suite
