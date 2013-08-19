'''
Tests various branching and search strategies. Compares number of LPs solved
for different branching and searching strategies.

The result of the script is as follows (until python decides to change its
random number generator).

problem     | branching strategy   | search strategy   | num lp
---------------------------------------------------------------
(10, 10, 0) | Most Fraction        | Depth First       | 7
(10, 10, 0) | Most Fraction        | Best First        | 5
(10, 10, 0) | Most Fraction        | Best Estimate     | 5
(10, 10, 0) | Fixed Branching      | Depth First       | 7
(10, 10, 0) | Fixed Branching      | Best First        | 5
(10, 10, 0) | Fixed Branching      | Best Estimate     | 5
(10, 10, 0) | Pseudocost Branching | Depth First       | 7
(10, 10, 0) | Pseudocost Branching | Best First        | 5
(10, 10, 0) | Pseudocost Branching | Best Estimate     | 5
(10, 10, 1) | Most Fraction        | Depth First       | 9
(10, 10, 1) | Most Fraction        | Best First        | 9
(10, 10, 1) | Most Fraction        | Best Estimate     | 9
(10, 10, 1) | Fixed Branching      | Depth First       | 9
(10, 10, 1) | Fixed Branching      | Best First        | 9
(10, 10, 1) | Fixed Branching      | Best Estimate     | 9
(10, 10, 1) | Pseudocost Branching | Depth First       | 9
(10, 10, 1) | Pseudocost Branching | Best First        | 9
(10, 10, 1) | Pseudocost Branching | Best Estimate     | 9
(20, 10, 2) | Most Fraction        | Depth First       | 37
(20, 10, 2) | Most Fraction        | Best First        | 33
(20, 10, 2) | Most Fraction        | Best Estimate     | 33
(20, 10, 2) | Fixed Branching      | Depth First       | 35
(20, 10, 2) | Fixed Branching      | Best First        | 31
(20, 10, 2) | Fixed Branching      | Best Estimate     | 31
(20, 10, 2) | Pseudocost Branching | Depth First       | 23
(20, 10, 2) | Pseudocost Branching | Best First        | 23
(20, 10, 2) | Pseudocost Branching | Best Estimate     | 23
(20, 10, 3) | Most Fraction        | Depth First       | 41
(20, 10, 3) | Most Fraction        | Best First        | 21
(20, 10, 3) | Most Fraction        | Best Estimate     | 24
(20, 10, 3) | Fixed Branching      | Depth First       | 57
(20, 10, 3) | Fixed Branching      | Best First        | 39
(20, 10, 3) | Fixed Branching      | Best Estimate     | 39
(20, 10, 3) | Pseudocost Branching | Depth First       | 31
(20, 10, 3) | Pseudocost Branching | Best First        | 21
(20, 10, 3) | Pseudocost Branching | Best Estimate     | 25
(30, 20, 4) | Most Fraction        | Depth First       | 95
(30, 20, 4) | Most Fraction        | Best First        | 95
(30, 20, 4) | Most Fraction        | Best Estimate     | 95
(30, 20, 4) | Fixed Branching      | Depth First       | 145
(30, 20, 4) | Fixed Branching      | Best First        | 145
(30, 20, 4) | Fixed Branching      | Best Estimate     | 145
(30, 20, 4) | Pseudocost Branching | Depth First       | 69
(30, 20, 4) | Pseudocost Branching | Best First        | 77
(30, 20, 4) | Pseudocost Branching | Best Estimate     | 83
(30, 20, 5) | Most Fraction        | Depth First       | 181
(30, 20, 5) | Most Fraction        | Best First        | 181
(30, 20, 5) | Most Fraction        | Best Estimate     | 181
(30, 20, 5) | Fixed Branching      | Depth First       | 159
(30, 20, 5) | Fixed Branching      | Best First        | 159
(30, 20, 5) | Fixed Branching      | Best Estimate     | 159
(30, 20, 5) | Pseudocost Branching | Depth First       | 113
(30, 20, 5) | Pseudocost Branching | Best First        | 111
(30, 20, 5) | Pseudocost Branching | Best Estimate     | 113
(40, 20, 6) | Most Fraction        | Depth First       | 323
(40, 20, 6) | Most Fraction        | Best First        | 181
(40, 20, 6) | Most Fraction        | Best Estimate     | 209
(40, 20, 6) | Fixed Branching      | Depth First       | 949
(40, 20, 6) | Fixed Branching      | Best First        | 525
(40, 20, 6) | Fixed Branching      | Best Estimate     | 562
(40, 20, 6) | Pseudocost Branching | Depth First       | 177
(40, 20, 6) | Pseudocost Branching | Best First        | 127
(40, 20, 6) | Pseudocost Branching | Best Estimate     | 151
(40, 20, 7) | Most Fraction        | Depth First       | 145
(40, 20, 7) | Most Fraction        | Best First        | 71
(40, 20, 7) | Most Fraction        | Best Estimate     | 74
(40, 20, 7) | Fixed Branching      | Depth First       | 111
(40, 20, 7) | Fixed Branching      | Best First        | 81
(40, 20, 7) | Fixed Branching      | Best Estimate     | 81
(40, 20, 7) | Pseudocost Branching | Depth First       | 81
(40, 20, 7) | Pseudocost Branching | Best First        | 55
(40, 20, 7) | Pseudocost Branching | Best Estimate     | 51
(40, 30, 8) | Most Fraction        | Depth First       | 691
(40, 30, 8) | Most Fraction        | Best First        | 411
(40, 30, 8) | Most Fraction        | Best Estimate     | 424
(40, 30, 8) | Fixed Branching      | Depth First       | 1083
(40, 30, 8) | Fixed Branching      | Best First        | 393
(40, 30, 8) | Fixed Branching      | Best Estimate     | 393
(40, 30, 8) | Pseudocost Branching | Depth First       | 395
(40, 30, 8) | Pseudocost Branching | Best First        | 167
(40, 30, 8) | Pseudocost Branching | Best Estimate     | 253

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


