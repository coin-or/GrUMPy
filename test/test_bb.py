'''
Tests correctness of branch and bound algorithm. Optimal values of problems are
hard-coded. This will be valid until python decides to change its random number
generator. See test_bb_pulp.py for comparing branch and bound results with
optimal values given by PuLP.

Script raises exceptions if the bb solution is not integer feasible or optimal
value is not right.
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
           (40,30,8),
           ]
# optimal values of problems
pre_computed_opt_val= {(10,10,0):50,
           (10,10,1):45,
           (20,10,2):104,
           (20,10,3):93,
           (30,20,4):139,
           (30,20,5):136,
           (40,20,6):192,
           (40,20,7):164,
           (40,30,8):179,
                       }

if __name__=='__main__':
    for p in problem:
        bt = BBTree()
        var, con, seed = p
        CONSTRAINTS, VARIABLES, OBJ, MAT, RHS = bt.GenerateRandomMIP(numVars=var,
                                                                     numCons=con,
                                                                     rand_seed=seed)
        solution, opt_value = bt.BranchAndBound(CONSTRAINTS, VARIABLES, OBJ,
                                                MAT, RHS)
        # test solution.
        #= test integer feasibility
        for v in solution:
            diff = solution[v]-math.floor(solution[v])
            if (diff>EPSILON and diff<(1-EPSILON)):
                raise Exception('Integer infeasible variable %s, value %f ' %(v, solution[v]))
        #= test feasibility of constraints
        Ax = []
        num_cons = len(CONSTRAINTS)
        #== for each constraint
        for c in range(num_cons):
            _sum = 0.0
            for v in VARIABLES:
                _sum += MAT[v][c]*solution[v]
            Ax.append(_sum)
        for c in range(num_cons):
            if Ax[c] > RHS[c]:
                raise Exception('Solution does not satisfy constraint ' + CONSTRAINTS[c])
        #= test optimal value
        if opt_value!=pre_computed_opt_val[p]:
            raise Exception('Optimality is not acheived for problem %s. BB: %f, OPT: %f ' %(str(p), opt_value[p], pre_computed_opt_val[p]))
        print '***************************************************'
        print '* No exceptions raised, BB solutions are correct. *'
        print '***************************************************'
