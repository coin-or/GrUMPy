from baktree import BAKTree
import pdb
import sys

if __name__ == '__main__':
#    pdb.set_trace()
    bt = BAKTree()
    bt.set_display_mode('pygame')
    while True:
        line = sys.stdin.readline()
        bt.ProcessLine(line)
        if bt.root is not None:
            bt.display()

