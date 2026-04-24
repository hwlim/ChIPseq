def check_column_counts(path, sep='\t'):
    """Check that all non-comment, non-empty lines have the same number of columns.
    Prints errors with line numbers. Returns True if file is consistent, False otherwise.
    """
    expected = None
    errors = []
    
    with open(path) as f:
        for lineno, line in enumerate(f, start=1):
            stripped = line.strip()
            # Skip empty lines and comment lines
            if not stripped or stripped.startswith('#'):
                continue
            
            n_cols = len(line.rstrip('\n').split(sep))
            
            if expected is None:
                expected = n_cols  # first data line sets the expectation
            elif n_cols != expected:
                errors.append((lineno, n_cols))
    
    if errors:
        print(f"ERROR: Expected {expected} columns (from first data line).")
        for lineno, n in errors:
            print(f"  Line {lineno}: {n} columns")
        return False
    else:
        print(f"OK: All data lines have {expected} columns.")
        return True


if not check_column_counts(src_sampleInfo):
    raise ValueError("TSV file has inconsistent column counts")
