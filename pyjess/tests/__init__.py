from . import test_tess_atom, test_test_template


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_tess_atom))
    suite.addTests(loader.loadTestsFromModule(test_test_template))
    # test_doctest.load_tests(loader, suite, pattern)
    return suite
