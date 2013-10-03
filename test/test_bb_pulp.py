'''
Tests correctness of branch and bound algorithm. Optimal values of problems are
compared to PuLP (a python interface/modeler to various open-source solvers)
results. See test_bb.py for comparing branch and bound results to hard-coded
values (no PuLP necessary).

The result of the script is as follows (until python decides to change its
random number generator).

Problem     | BB optimal   | PuLP optimal
-----------------------------------------
(10, 10, 0) | 50.0         | 50.0
(10, 10, 1) | 45.0         | 45.0
(20, 10, 2) | 104.0        | 104.0
(20, 10, 3) | 93.0         | 93.0
(30, 20, 4) | 139.0        | 139.0
(30, 20, 5) | 136.0        | 136.0
(40, 20, 6) | 192.0        | 192.0
(40, 20, 7) | 164.0        | 164.0
(40, 30, 8) | 179.0        | 179.0

'''

from grumpy import BBTree
# import branching strategies
from grumpy import MOST_FRACTIONAL, FIXED_BRANCHING, PSEUDOCOST_BRANCHING
# import searching strategies
from grumpy import DEPTH_FIRST, BEST_FIRST, BEST_ESTIMATE
import pulp
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
# will keep PuLP optimal objective value
pulp_optimal = {}
# will keep branch and bound optimal objective value
bb_optimal = {}

if __name__=='__main__':
    for p in problem:
        bt = BBTree()
        var, con, seed = p
        CONSTRAINTS, VARIABLES, OBJ, MAT, RHS = bt.GenerateRandomMIP(numVars=var,
                                                                     numCons=con,
                                                                     rand_seed=seed)
        # construct pulp problem
        # create problem object
        prob = pulp.LpProblem('test', sense=pulp.LpMaximize)
        # create variable dictionary
        prob_var = {}
        for v in VARIABLES:
            prob_var[v] = pulp.LpVariable(v, cat = pulp.LpBinary)
        # add constraints
        for c in range(con):
            prob += pulp.lpSum([MAT[v][c]*prob_var[v] for v in VARIABLES]) <= RHS[c], CONSTRAINTS[c]
        # add objective
        obj = pulp.lpSum([OBJ[v]*prob_var[v] for v in VARIABLES])
        prob += obj, 'objective'
        # solve pulp problem
        status = prob.solve()
        # check status
        if pulp.LpStatus[status] != 'Optimal':
            raise Exception('Pulp: Optimal solution could not be found.')
        # get optimal value
        pulp_optimal[p] = pulp.value(obj)
        # solve using BB
        solution, bb_optimal[p] = bt.BranchAndBound(CONSTRAINTS, VARIABLES, OBJ,
                                                 MAT, RHS)
        # test solution.
        #= test integer feasibility
        for v in solution:
            diff = solution[v]-math.floor(solution[v])
            if (diff>EPSILON and diff<(1-EPSILON)):
                print diff
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
        # test optimality
        if bb_optimal[p] != pulp_optimal[p]:
            raise Exception('Optimality is not acheived for problem %s. BB: %f, OPT: %f ' %(str(p), bb_optimal, pulp_optimal))
    # print optimal values in a table
    print 'Problem     | BB optimal   | PuLP optimal'
    print '-----------------------------------------'
    for p in problem:
        print str(p).ljust(10),
        print '|',
        print str(bb_optimal[p]).ljust(12),
        print '|',
        print str(pulp_optimal[p])
