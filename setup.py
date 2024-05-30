import configparser
import functools
import glob
import io
import multiprocessing.pool
import os
import platform
import re
import setuptools
import sys
import sysconfig
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Utils ------------------------------------------------------------------


def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)


def _patch_osx_compiler(compiler):
    # On newer OSX, Python has been compiled as a universal binary, so
    # it will attempt to pass universal binary flags when building the
    # extension. This will not work because the code makes use of SSE2.
    for tool in ("compiler", "compiler_so", "linker_so"):
        flags = getattr(compiler, tool)
        i = next(
            (
                i
                for i in range(1, len(flags))
                if flags[i - 1] == "-arch" and flags[i] != platform.machine()
            ),
            None,
        )
        if i is not None:
            flags.pop(i)
            flags.pop(i - 1)


def _detect_target_machine(platform):
    if platform == "win32":
        return "x86"
    return platform.rsplit("-", 1)[-1]


def _detect_target_cpu(platform):
    machine = _detect_target_machine(platform)
    if re.match("^mips", machine):
        return "mips"
    elif re.match("^(aarch64|arm64)$", machine):
        return "aarch64"
    elif re.match("^arm", machine):
        return "arm"
    elif re.match("(x86_64)|AMD64|amd64", machine):
        return "x86_64"
    elif re.match("(x86)|(^i.86$)", machine):
        return "x86"
    elif re.match("^(powerpc|ppc)", machine):
        return "ppc"
    return None


def _detect_target_system(platform):
    if platform.startswith("win"):
        return "windows"
    elif platform.startswith("macos"):
        return "macos"
    elif platform.startswith("linux"):
        return "linux_or_android"
    elif platform.startswith("freebsd"):
        return "freebsd"
    return None


# --- Commands ------------------------------------------------------------------


class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly."""

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", "build-backend", '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class build_clib(_build_clib):
    """A custom `build_clib` that makes all C++ class attributes public."""

    # --- Compatibility with `setuptools.Command`

    user_options = _build_clib.user_options + [
        ("parallel", "j", "number of parallel build jobs"),
    ]

    def initialize_options(self):
        _build_clib.initialize_options(self)
        self.parallel = None

    def finalize_options(self):
        _build_clib.finalize_options(self)
        if self.parallel is not None:
            self.parallel = int(self.parallel)

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [lib.name for lib in self.libraries]

    def get_source_files(self):
        return [source for lib in self.libraries for source in lib.sources]

    def get_library(self, name):
        return next(lib for lib in self.libraries if lib.name == name)

    # --- Patch and copy sources ---

    _HEADER_PATTERN = re.compile(r"^@@ -(\d+),?(\d+)? \+(\d+),?(\d+)? @@.*$")

    @classmethod
    def _apply_patch(cls, s, patch, revert=False):
        # see https://stackoverflow.com/a/40967337
        s = s.splitlines(keepends=True)
        p = patch.splitlines(keepends=True)[2:]  # ignore `git diff` header
        t = []
        i = 0
        sl = 0
        midx, sign = (1, "+") if not revert else (3, "-")
        while i < len(p) and p[i].startswith(("---", "+++")):
            i += 1  # skip header lines

        while i < len(p):
            match = cls._HEADER_PATTERN.match(p[i])
            if not match:
                raise ValueError("Invalid line in patch: {!r}".format(p[i]))
            i += 1
            l = int(match.group(midx)) - 1 + (match.group(midx + 1) == "0")
            t.extend(s[sl:l])
            sl = l
            while i < len(p) and p[i][0] != "@":
                if i + 1 < len(p) and p[i + 1][0] == "\\":
                    line = p[i][:-1]
                    i += 2
                else:
                    line = p[i]
                    i += 1
                if len(line) > 0:
                    if line[0] == sign or line[0] == " ":
                        t += line[1:]
                    sl += line[0] != sign

        t.extend(s[sl:])
        return "".join(t)

    def _patch_file(self, input, output):
        basename = os.path.basename(input)
        patchname = os.path.realpath(
            os.path.join(__file__, os.pardir, "patches", "{}.patch".format(basename))
        )
        if os.path.exists(patchname):
            _eprint(
                "patching", os.path.relpath(input), "with", os.path.relpath(patchname)
            )
            with open(patchname, "r") as patchfile:
                patch = patchfile.read()
            with open(input, "r") as src:
                srcdata = src.read()
            with open(output, "w") as dst:
                dst.write(self._apply_patch(srcdata, patch))
        else:
            self.copy_file(input, output)

    # --- Build code ---

    def build_libraries(self, libraries):
        # remove universal compilation flags for OSX
        if platform.system() == "Darwin":
            _patch_osx_compiler(self.compiler)

        # build each library only if the sources are outdated
        self.mkpath(self.build_clib)
        for library in libraries:
            libname = self.compiler.library_filename(
                library.name, output_dir=self.build_clib
            )
            self.make_file(library.sources, libname, self.build_library, (library,))

    def build_library(self, library):
        # show the compiler being used
        _eprint(
            "building", library.name, "with", self.compiler.compiler_type, "compiler"
        )

        # add debug flags if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-g")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/Z7")

        # add Windows flags
        if self.compiler.compiler_type == "msvc":
            library.define_macros.append(("WIN32", 1))

        # add Linux flags
        if platform.system() == "Linux":
            library.define_macros.append(("_GNU_SOURCE", 1))

        # copy and patch headers to build directory
        for header in library.depends:
            output = os.path.join(self.build_clib, os.path.basename(header))
            self.mkpath(os.path.dirname(output))
            self.make_file([header], output, self._patch_file, (header, output))

        # copy and patch sources to build directory
        sources = [
            os.path.join(self.build_temp, os.path.basename(source))
            for source in library.sources
        ]
        for source, source_copy in zip(library.sources, sources):
            self.make_file(
                [source], source_copy, self._patch_file, (source, source_copy)
            )

        # store compile args
        compile_args = (
            None,
            library.define_macros,
            library.include_dirs + [self.build_clib],
            self.debug,
            library.extra_compile_args,
            None,
            library.depends,
        )
        # manually prepare sources and get the names of object files
        objects = [
            re.sub(r"(.cpp|.c)$", self.compiler.obj_extension, s) for s in sources
        ]
        # only compile outdated files
        with multiprocessing.pool.ThreadPool(self.parallel) as pool:
            pool.starmap(
                functools.partial(self._compile_file, compile_args=compile_args),
                zip(sources, objects),
            )

        # link into a static library
        libfile = self.compiler.library_filename(
            library.name,
            output_dir=self.build_clib,
        )
        self.make_file(
            objects,
            libfile,
            self.compiler.create_static_lib,
            (objects, library.name, self.build_clib, None, self.debug),
        )

    def _compile_file(self, source, object, compile_args):
        self.make_file(
            [source], object, self.compiler.compile, ([source], *compile_args)
        )


class build_ext(_build_ext):
    """A `build_ext` that adds various flags and defines."""

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.target_machine = None
        self.target_system = None
        self.target_cpu = None

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # check platform
        if self.plat_name is None:
            self.plat_name = sysconfig.get_platform()
        # detect platform options
        self.target_machine = _detect_target_machine(self.plat_name)
        self.target_system = _detect_target_system(self.plat_name)
        self.target_cpu = _detect_target_cpu(self.plat_name)
        # transfer arguments to the build_clib method
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.debug = self.debug
        self._clib_cmd.force = self.force
        self._clib_cmd.verbose = self.verbose
        self._clib_cmd.define = self.define
        self._clib_cmd.include_dirs = self.include_dirs
        self._clib_cmd.compiler = self.compiler
        self._clib_cmd.parallel = self.parallel

    def _check_getid(self):
        _eprint("checking whether `PyInterpreterState_GetID` is available")

        base = "have_getid"
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        objects = []

        self.mkpath(self.build_temp)
        with open(testfile, "w") as f:
            f.write(
                """
            #include <stdint.h>
            #include <stdlib.h>
            #include <Python.h>

            int main(int argc, char *argv[]) {{
                PyInterpreterState_GetID(NULL);
                return 0;
            }}
            """
            )

        if self.compiler.compiler_type == "msvc":
            flags = ["/WX"]
        else:
            flags = ["-Werror=implicit-function-declaration"]

        try:
            self.mkpath(self.build_temp)
            objects = self.compiler.compile([testfile], extra_postargs=flags)
        except CompileError:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)

    def build_extension(self, ext):
        # show the compiler being used
        _eprint(
            "building",
            ext.name,
            "with",
            self.compiler.compiler_type,
            "compiler for platform",
            self.plat_name,
        )

        # add debug symbols if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-g")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Z7")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        else:
            ext.define_macros.append(("CYTHON_WITHOUT_ASSERTIONS", 1))

        # check if `PyInterpreterState_GetID` is defined
        if self._check_getid():
            ext.define_macros.append(("HAS_PYINTERPRETERSTATE_GETID", 1))

        # use the compiled library
        ext.include_dirs.append(self._clib_cmd.build_clib)
        for libname in ext.libraries:
            library = self._clib_cmd.get_library(libname)
            ext.extra_objects.append(
                self.compiler.library_filename(
                    library.name,
                    output_dir=self._clib_cmd.build_clib,
                )
            )

        # build the rest of the extension as normal
        ext._needs_stub = False
        _build_ext.build_extension(self, ext)

    def build_extensions(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError(
                "Cython is required to run `build_ext` command"
            ) from cythonize

        # check build_clib has been run
        if not self.distribution.have_run.get("build_clib", False):
            self._clib_cmd.run()

        # remove universal compilation flags for OSX
        if platform.system() == "Darwin":
            _patch_osx_compiler(self.compiler)

        # use debug directives with Cython if building in debug mode
        cython_args = {
            "include_path": ["include"],
            "compiler_directives": {
                "cdivision": True,
                "nonecheck": False,
            },
            "compile_time_env": {
                "SYS_IMPLEMENTATION_NAME": sys.implementation.name,
                "SYS_VERSION_INFO_MAJOR": sys.version_info.major,
                "SYS_VERSION_INFO_MINOR": sys.version_info.minor,
                "SYS_VERSION_INFO_MICRO": sys.version_info.micro,
                "MAX_ALPHABET_SIZE": 32,
                "DEFAULT_BUFFER_SIZE": io.DEFAULT_BUFFER_SIZE,
                "TARGET_CPU": self.target_cpu,
                "TARGET_SYSTEM": self.target_system,
            },
        }
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["cdivision_warnings"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
        else:
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False

        # cythonize the extensions
        self.extensions = cythonize(self.extensions, **cython_args)

        # build the extensions as normal
        _build_ext.build_extensions(self)


class clean(_clean):
    """A `clean` that removes intermediate files created by Cython."""

    def run(self):

        source_dir = os.path.join(os.path.dirname(__file__), "pyswrd")

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                _eprint("removing {!r}".format(file))
                os.remove(file)

        _clean.run(self)


# --- Setup ---------------------------------------------------------------------

setuptools.setup(
    libraries=[
        Library(
            "jess",
            language="c",
            include_dirs=[os.path.join("vendor", "jess")],
            libraries=["m"],
            sources=[
                os.path.join("vendor", "jess", "src", "Annulus.c"),
                os.path.join("vendor", "jess", "src", "Atom.c"),
                os.path.join("vendor", "jess", "src", "Jess.c"),
                os.path.join("vendor", "jess", "src", "Join.c"),
                os.path.join("vendor", "jess", "src", "KdTree.c"),
                os.path.join("vendor", "jess", "src", "Molecule.c"),
                os.path.join("vendor", "jess", "src", "Region.c"),
                os.path.join("vendor", "jess", "src", "Scanner.c"),
                os.path.join("vendor", "jess", "src", "Super.c"),
                os.path.join("vendor", "jess", "src", "TessAtom.c"),
                os.path.join("vendor", "jess", "src", "TessTemplate.c"),
            ],
            depends=[
                os.path.join("vendor", "jess", "src", "Annulus.h"),
                os.path.join("vendor", "jess", "src", "Atom.h"),
                os.path.join("vendor", "jess", "src", "Jess.h"),
                os.path.join("vendor", "jess", "src", "Join.h"),
                os.path.join("vendor", "jess", "src", "KdTree.h"),
                os.path.join("vendor", "jess", "src", "Molecule.h"),
                os.path.join("vendor", "jess", "src", "Region.h"),
                os.path.join("vendor", "jess", "src", "Scanner.h"),
                os.path.join("vendor", "jess", "src", "Super.h"),
                os.path.join("vendor", "jess", "src", "Template.h"),
                os.path.join("vendor", "jess", "src", "TessAtom.h"),
                os.path.join("vendor", "jess", "src", "TessTemplate.h"),
            ],
        ),
    ],
    ext_modules=[
        Extension(
            "pyjess._jess",
            language="c",
            libraries=["jess"],
            sources=[
                os.path.join("pyjess", "_jess.pyx"),
            ],
            include_dirs=[
                "pyjess",
                "include",
            ],
        ),
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "build_clib": build_clib,
        "clean": clean,
    },
)
