'''
Tests correctness of branch and bound algorithm.
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
           (40,30,8)
           ]

if __name__=='__main__':
    for p in problem:
        bt = BBTree()
        var, con, seed = p
        CONSTRAINTS, VARIABLES, OBJ, MAT, RHS = bt.GenerateRandomMIP(numVars=var,
                                                                     numCons=con,
                                                                     rand_seed=seed)
        # construct pulp problem
        # create problem object
        prob = pulp.LpProblem('test', sense=pulp.LpMinimize)
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
        pulp_optimal = pulp.value(obj)
        # solve using BB
        solution, bb_optimal = bt.BranchAndBound(CONSTRAINTS, VARIABLES, OBJ,
                                                 MAT, RHS)
        # test solution.
        #= test integer feasibility
        for v in solution:
            diff = solution[v]-math.floor(solution[v])
            if (diff>EPSILON and diff<(1-EPSILON)):
                print diff
                raise Exception('Integer infeasible variable %s, value %f ' %(v, solution[v]))
        # test optimality
        if bb_optimal < pulp_optimal:
            raise Exception('Optimality is not acheived for problem %s. BB: %f, OPT: %f ' %(str(p), bb_optimal, pulp_optimal))

