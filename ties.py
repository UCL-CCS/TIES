#!/home/dresio/software/virtualenvs/ties/bin/python
# -*- coding: utf-8 -*-
"""
This is a helpful file for the TIES interface for debugging in Pycharm. It can be run without installation and
relays the commands to the ties/ties.py, the same way it is done when it is installed by pip in the python package.
"""
import re
import sys
from ties import ties
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(ties.command_line_script())
(ties)
