'''
How to control display sequence on xdot:
1. Give a valid display_interval. Display interval represents the number of
nodes that will be added before next display call, ie. next tree display on
your gtk window. When display method is called BB stops. It will resume when
you close the gtk window. Gtk window will pop up again when the required number
of nodes (display_interval many) is added to the tree.
2. If you want to get only the final tree give sys.maxiter as input.

Examples for display modes/layouts:
1. file display mode, no layout specified, save in png format
bt = BBTree(display='file')

2. xdot display mode
bt = BBTree(display='xdot')

3. pygame display mode, no layout specified, save in png format
bt = BBTree(display='pygame')

4. file display mode, with dot2tex layout
bt = BBTree(display='file', layout = 'dot2tex')

'''
from grumpy import BBTree
import sys

if __name__=='__main__':
    #======== file display mode, no layout specified, save in png format
    #bt = BBTree(display='file')
    #======== xdot display mode
    #bt = BBTree(display='xdot')
    #======== pygame display mode, no layout specified, save in png format
    bt = BBTree(display='pygame')
    #======== file display mode, with dot2tex layout
    #bt = BBTree(display='file', layout = 'dot2tex')

    CONSTRAINTS, VARIABLES, OBJ, MAT, RHS = bt.GenerateRandomMIP(numVars=30, numCons=10, rand_seed = 0)
    bt.BranchAndBound(CONSTRAINTS, VARIABLES, OBJ, MAT, RHS,
                     branch_strategy = 'Pseudocost', search_strategy = 'Best First',
                     display_interval = 4)

