#!/home/dresio/software/virtualenvs/ties/bin/python
# -*- coding: utf-8 -*-
import re
import sys
from ties import ties
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(ties.command_line_script())
(ties)
