'''
Tests various branching and search strategies.
'''

from grumpy import BBTree
# import branching strategies
from grumpy import MOST_FRACTIONAL, FIXED_BRANCHING, PSEUDOCOST_BRANCHING
# import searching strategies
from grumpy import DEPTH_FIRST, BEST_FIRST, BEST_ESTIMATE
import sys
import math

branch_strategy = [MOST_FRACTIONAL, FIXED_BRANCHING, PSEUDOCOST_BRANCHING]
search_strategy = [DEPTH_FIRST, BEST_FIRST, BEST_ESTIMATE]
# test problem, (num_vars,num_cons,seed)
problem = [(10,10,0),
           (10,10,1),
           (20,10,2),
           (20,10,3),
           (30,20,4),
           (30,20,5),
           (40,20,6),
           (40,20,7),
           (40,30,8)
           ]
# number of LPs solved for each problem
num_lp = {}

if __name__=='__main__':
    for p in problem:
        var, con, seed = p
        bt = BBTree()
        CONSTRAINTS, VARIABLES, OBJ, MAT, RHS = bt.GenerateRandomMIP(numVars=var,
                                                                     numCons=con,
                                                                     rand_seed=seed)
        num_lp[p] = {}
        for b in branch_strategy:
            for s in search_strategy:
                # create a new object before each branch and bound call
                bt = BBTree()
                # solve using BB
                solution, bb_optimal = bt.BranchAndBound(CONSTRAINTS, VARIABLES, OBJ,
                                                         MAT, RHS,
                                                         branch_strategy = b,
                                                         search_strategy = s)
                # keep number of LPs solved.
                num_lp[p][(b,s)] = bt._lp_count
    # print table
    print 'problem     | branching strategy   | search strategy   | num lp'
    print '---------------------------------------------------------------'
    for p in problem:
        for b in branch_strategy:
            for s in search_strategy:
                print str(p).ljust(10),
                print '|',
                print str(b).ljust(20),
                print '|',
                print str(s).ljust(17),
                print '|',
                print num_lp[p][(b,s)]


