#!/usr/bin/env python

import pytest, tempfile, os, sys, site, shutil

def main():
    with tempfile.TemporaryDirectory() as tmpdirname:
        print('created temporary directory', tmpdirname)
        test_dir = os.path.join(site.getsitepackages()[0], 'ties', 'testing')
        tmp_test_dir = os.path.join(tmpdirname, 'ties', 'unit_testing', 'testing')
        shutil.copytree(test_dir, tmp_test_dir)
        sys.exit(pytest.main([os.path.join(tmp_test_dir)]))

if __name__ == '__main__':
    main()
