try:
    from importlib.resources import files
except ImportError:
    try:
        from importlib_resources import files # type: ignore
    except ImportError:
        files = None # type: ignore