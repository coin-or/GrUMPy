'''
Tests correctness of branch and bound algorithm.
'''

from grumpy import BBTree
import sys
import math

EPSILON = 1e-10
# test problem, (num_vars,num_cons,seed)
problem = [(10,10,0),
           (10,10,1),
           (20,10,2),
           (20,10,3),
           (30,20,4),
           (30,20,5),
           (40,20,6),
           (40,20,7),
           (40,30,8)]
# optimal values of problems
pre_computed_opt_val= {(10,10,0):,
           (10,10,1):,
           (20,10,2):,
           (20,10,3):,
           (30,20,4):,
           (30,20,5):,
           (40,20,6):,
           (40,20,7):,
           (40,30,8): }

if __name__=='__main__':
    for p in problem:
        bt = BBTree()
        numVars, numCons, seed = p
        CONSTRAINTS, VARIABLES, OBJ, MAT, RHS = bt.GenerateRandomMIP(numVars=var,
                                                                     numCons=con,
                                                                     rand_seed=seed)
        solution, opt_value = bt.BranchAndBound(CONSTRAINTS, VARIABLES, OBJ,
                                                MAT, RHS,
                                                branch_strategy = 'Pseudocost',
                                                search_strategy = 'Best First')
        # test solution.
        #= test integer feasibility
        for v in solution:
            diff = v-math.floor(v)
            if (diff>EPSILON or diff<(1-EPSILON)):
                raise Exception('Integer infeasible variable %d, value %f ' %(i, solution[i]))
        if opt_value<pre_computed_opt_val[p]:
            raise Exception('Optimality is not acheived for problem %s. BB: %f, OPT: %f ' %(str(p), opt_value[p], pre_computed_opt_val[p]))

