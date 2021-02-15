from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from .BBTree import *
from .BranchAndBound import *
try:
    from .polyhedron2D import *
except ImportError:
    pass
import sys
import io

if __name__ == '__main__':
    bt = BBTree()
    bt.set_display_mode('matplotlib')
    line_number = 0
    for line in sys.stdin:
        bt.ProcessLine(line)
        print('line', line_number, 'processed.')
        line_number = line_number+1
        if line_number%10 != 0:
            continue
        bt.display_all()
