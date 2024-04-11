from . import test_atom, test_tess_atom, test_tess_template


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_atom))
    suite.addTests(loader.loadTestsFromModule(test_tess_atom))
    suite.addTests(loader.loadTestsFromModule(test_tess_template))
    # test_doctest.load_tests(loader, suite, pattern)
    return suite
