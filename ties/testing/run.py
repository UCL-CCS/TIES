#!/usr/bin/env python

import pytest, tempfile, os, sys, site, shutil
from pathlib import Path


def main():
    try:
        with tempfile.TemporaryDirectory() as tmpdirname:
            cwd = Path.cwd()
            os.chdir(tmpdirname)
            test_dir = os.path.join(site.getsitepackages()[0], "ties", "testing")
            examples = os.path.join(site.getsitepackages()[0], "ties", "examples")
            shutil.copytree(test_dir, "testing")
            shutil.copytree(examples, "examples")

            # enter the testing directory
            os.chdir(Path(tmpdirname) / "testing")
            sys.exit(pytest.main([Path(tmpdirname) / "testing"]))
    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    main()
