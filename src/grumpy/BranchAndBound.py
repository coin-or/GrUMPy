from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from past.utils import old_div
__author__ = 'Ted Ralphs'
__maintainer__ = 'Ted Ralphs (ted@lehigh.edu)'

import random, sys, math
try:
    from src.blimpy import PriorityQueue
except ImportError:
    from coinor.blimpy import PriorityQueue
import time
from pulp import LpVariable, lpSum, LpProblem, LpMaximize, LpConstraint
from pulp import LpStatus, value
from .BBTree import BBTree
from .BBTree import MOST_FRACTIONAL, FIXED_BRANCHING, PSEUDOCOST_BRANCHING
from .BBTree import DEPTH_FIRST, BEST_FIRST, BEST_ESTIMATE, INFINITY

def GenerateRandomMIP(numVars = 40, numCons = 20, density = 0.2,
                      maxObjCoeff = 10, maxConsCoeff = 10, 
                      tightness = 2, rand_seed = 2, layout = 'dot'):
    random.seed(rand_seed)
    CONSTRAINTS = ["C"+str(i) for i in range(numCons)]
    if layout == 'dot2tex':
        VARIABLES = ["x_{"+str(i)+"}" for i in range(numVars)]
    else:
        VARIABLES = ["x"+str(i) for i in range(numVars)]
    OBJ = dict((i, random.randint(1, maxObjCoeff)) for i in VARIABLES)
    MAT = dict((i, [random.randint(1, maxConsCoeff)
                    if random.random() <= density else 0
                    for j in CONSTRAINTS]) for i in VARIABLES)
    RHS = [random.randint(int(numVars*density*maxConsCoeff/tightness),
                   int(numVars*density*maxConsCoeff/1.5))
           for i in CONSTRAINTS]
    return CONSTRAINTS, VARIABLES, OBJ, MAT, RHS

def BranchAndBound(T, CONSTRAINTS, VARIABLES, OBJ, MAT, RHS,
                   branch_strategy = MOST_FRACTIONAL,
                   search_strategy = DEPTH_FIRST,
                   complete_enumeration = False,
                   display_interval = None,
                   binary_vars = True):
    
    if T.get_layout() == 'dot2tex':
        cluster_attrs = {'name':'Key', 'label':r'\text{Key}', 'fontsize':'12'}
        T.add_node('C', label = r'\text{Candidate}', style = 'filled',
                      color = 'yellow', fillcolor = 'yellow')
        T.add_node('I', label = r'\text{Infeasible}', style = 'filled',
                      color = 'orange', fillcolor = 'orange')
        T.add_node('S', label = r'\text{Solution}', style = 'filled',
                      color = 'lightblue', fillcolor = 'lightblue')
        T.add_node('P', label = r'\text{Pruned}', style = 'filled',
                      color = 'red', fillcolor = 'red')
        T.add_node('PC', label = r'\text{Pruned}$\\ $\text{Candidate}', style = 'filled',
                      color = 'red', fillcolor = 'yellow')
    else:
        cluster_attrs = {'name':'Key', 'label':'Key', 'fontsize':'12'}
        T.add_node('C', label = 'Candidate', style = 'filled',
                      color = 'yellow', fillcolor = 'yellow')
        T.add_node('I', label = 'Infeasible', style = 'filled',
                      color = 'orange', fillcolor = 'orange')
        T.add_node('S', label = 'Solution', style = 'filled',
                      color = 'lightblue', fillcolor = 'lightblue')
        T.add_node('P', label = 'Pruned', style = 'filled',
                      color = 'red', fillcolor = 'red')
        T.add_node('PC', label = 'Pruned \n Candidate', style = 'filled',
                      color = 'red', fillcolor = 'yellow')
    T.add_edge('C', 'I', style = 'invisible', arrowhead = 'none')
    T.add_edge('I', 'S', style = 'invisible', arrowhead = 'none')
    T.add_edge('S', 'P', style = 'invisible', arrowhead = 'none')
    T.add_edge('P', 'PC', style = 'invisible', arrowhead = 'none')
    T.create_cluster(['C', 'I', 'S', 'P', 'PC'], cluster_attrs)
    # The initial lower bound
    LB = -INFINITY
    # The number of LP's solved, and the number of nodes solved
    node_count = 1
    iter_count = 0
    lp_count = 0
    
    if binary_vars:
        var   = LpVariable.dicts("", VARIABLES, 0, 1)
    else:
        var   = LpVariable.dicts("", VARIABLES)
    
    numCons = len(CONSTRAINTS)
    numVars = len(VARIABLES)
    # List of incumbent solution variable values
    opt = dict([(i, 0) for i in VARIABLES])
    pseudo_u = dict((i, (OBJ[i], 0)) for i in VARIABLES)
    pseudo_d = dict((i, (OBJ[i], 0)) for i in VARIABLES)
    print("===========================================")
    print("Starting Branch and Bound")
    if branch_strategy == MOST_FRACTIONAL:
        print("Most fractional variable")
    elif branch_strategy == FIXED_BRANCHING:
        print("Fixed order")
    elif branch_strategy == PSEUDOCOST_BRANCHING:
        print("Pseudocost brancing")
    else:
        print("Unknown branching strategy %s" %branch_strategy)
    if search_strategy == DEPTH_FIRST:
        print("Depth first search strategy")
    elif search_strategy == BEST_FIRST:
        print("Best first search strategy")
    else:
        print("Unknown search strategy %s" %search_strategy)
    print("===========================================")
    # List of candidate nodes
    Q = PriorityQueue()
    # The current tree depth
    cur_depth = 0
    cur_index = 0
    # Timer
    timer = time.time()
    Q.push(0, -INFINITY, (0, None, None, None, None, None, None))
    # Branch and Bound Loop
    while not Q.isEmpty():
        infeasible = False
        integer_solution = False
        (cur_index, parent, relax, branch_var, branch_var_value, sense,
        rhs) = Q.pop()
        if cur_index is not 0:
            cur_depth = T.get_node_attr(parent, 'level') + 1
        else:
            cur_depth = 0
        print("")
        print("----------------------------------------------------")
        print("")
        if LB > -INFINITY:
            print("Node: %s, Depth: %s, LB: %s" %(cur_index,cur_depth,LB))
        else:
            print("Node: %s, Depth: %s, LB: %s" %(cur_index,cur_depth,"None"))
        if relax is not None and relax <= LB:
            print("Node pruned immediately by bound")
            T.set_node_attr(parent, 'color', 'red')
            continue
        #====================================
        #    LP Relaxation
        #====================================
        # Compute lower bound by LP relaxation
        prob = LpProblem("relax", LpMaximize)
        prob += lpSum([OBJ[i]*var[i] for i in VARIABLES]), "Objective"
        for j in range(numCons):
            prob += (lpSum([MAT[i][j]*var[i] for i in VARIABLES])<=RHS[j],\
                         CONSTRAINTS[j])
        # Fix all prescribed variables
        branch_vars = []
        if cur_index is not 0:
            sys.stdout.write("Branching variables: ")
            branch_vars.append(branch_var)
            if sense == '>=':
                prob += LpConstraint(lpSum(var[branch_var]) >= rhs)
            else:
                prob += LpConstraint(lpSum(var[branch_var]) <= rhs)
            print(branch_var, end=' ')
            pred = parent
            while not str(pred) == '0':
                pred_branch_var = T.get_node_attr(pred, 'branch_var')
                pred_rhs = T.get_node_attr(pred, 'rhs')
                pred_sense = T.get_node_attr(pred, 'sense')
                if pred_sense == '<=':
                    prob += LpConstraint(lpSum(var[pred_branch_var])
                                         <= pred_rhs)
                else:
                    prob += LpConstraint(lpSum(var[pred_branch_var])
                                         >= pred_rhs)
                print(pred_branch_var, end=' ')
                branch_vars.append(pred_branch_var)
                pred = T.get_node_attr(pred, 'parent')
            print()
        # Solve the LP relaxation
        prob.solve()
        lp_count = lp_count +1
        # Check infeasibility
        infeasible = LpStatus[prob.status] == "Infeasible" or \
            LpStatus[prob.status] == "Undefined"
        # Print status
        if infeasible:
            print("LP Solved, status: Infeasible")
        else:
            print("LP Solved, status: %s, obj: %s" %(LpStatus[prob.status],
                                                     value(prob.objective)))
        if(LpStatus[prob.status] == "Optimal"):
            relax = value(prob.objective)
            # Update pseudocost
            if branch_var != None:
                if sense == '<=':
                    pseudo_d[branch_var] = (
                    old_div((pseudo_d[branch_var][0]*pseudo_d[branch_var][1] +
                    old_div((T.get_node_attr(parent, 'obj') - relax),
                    (branch_var_value - rhs))),(pseudo_d[branch_var][1]+1)),
                    pseudo_d[branch_var][1]+1)
                else:
                    pseudo_u[branch_var] = (
                    old_div((pseudo_u[branch_var][0]*pseudo_d[branch_var][1] +
                     old_div((T.get_node_attr(parent, 'obj') - relax),
                     (rhs - branch_var_value))),(pseudo_u[branch_var][1]+1)),
                    pseudo_u[branch_var][1]+1)
            var_values = dict([(i, var[i].varValue) for i in VARIABLES])
            integer_solution = 1
            for i in VARIABLES:
                if (abs(round(var_values[i]) - var_values[i]) > .001):
                    integer_solution = 0
                    break
            # Determine integer_infeasibility_count and
            # Integer_infeasibility_sum for scatterplot and such
            integer_infeasibility_count = 0
            integer_infeasibility_sum = 0.0
            for i in VARIABLES:
                if (var_values[i] not in set([0,1])):
                    integer_infeasibility_count += 1
                    integer_infeasibility_sum += min([var_values[i],
                                                      1.0-var_values[i]])
            if (integer_solution and relax>LB):
                LB = relax
                for i in VARIABLES:
                    # These two have different data structures first one
                    #list, second one dictionary
                    opt[i] = var_values[i]
                print("New best solution found, objective: %s" %relax)
                for i in VARIABLES:
                    if var_values[i] > 0:
                        print("%s = %s" %(i, var_values[i]))
            elif (integer_solution and relax<=LB):
                print("New integer solution found, objective: %s" %relax)
                for i in VARIABLES:
                    if var_values[i] > 0:
                        print("%s = %s" %(i, var_values[i]))
            else:
                print("Fractional solution:")
                for i in VARIABLES:
                    if var_values[i] > 0:
                        print("%s = %s" %(i, var_values[i]))
            #For complete enumeration
            if complete_enumeration:
                relax = LB - 1
        else:
            relax = INFINITY
        if integer_solution:
            print("Integer solution")
            BBstatus = 'S'
            status = 'integer'
            color = 'lightblue'
        elif infeasible:
            print("Infeasible node")
            BBstatus = 'I'
            status = 'infeasible'
            color = 'orange'
        elif not complete_enumeration and relax <= LB:
            print("Node pruned by bound (obj: %s, UB: %s)" %(relax,LB))
            BBstatus = 'P'
            status = 'fathomed'
            color = 'red'
        elif cur_depth >= numVars :
            print("Reached a leaf")
            BBstatus = 'fathomed'
            status = 'L'
        else:
            BBstatus = 'C'
            status = 'candidate'
            color = 'yellow'
        if BBstatus is 'I':
            if T.get_layout() == 'dot2tex':
                label = '\text{I}'
            else:
                label = 'I'
        else:
            label = "%.1f"%relax
        if iter_count == 0:
            if status is not 'candidate':
                integer_infeasibility_count = None
                integer_infeasibility_sum = None
            if status is 'fathomed':
                if T._incumbent_value is None:
                    print('WARNING: Encountered "fathom" line before '+\
                        'first incumbent.')
            T.AddOrUpdateNode(0, None, None, 'candidate', relax,
                             integer_infeasibility_count,
                             integer_infeasibility_sum,
                             label = label,
                             obj = relax, color = color,
                             style = 'filled', fillcolor = color)
            if status is 'integer':
                T._previous_incumbent_value = T._incumbent_value
                T._incumbent_value = relax
                T._incumbent_parent = -1
                T._new_integer_solution = True
#           #Currently broken
#           if ETREE_INSTALLED and T.attr['display'] == 'svg':
#               T.write_as_svg(filename = "node%d" % iter_count,
#                                 nextfile = "node%d" % (iter_count + 1),
#                                 highlight = cur_index)
        else:
            _direction = {'<=':'L', '>=':'R'}
            if status is 'infeasible':
                integer_infeasibility_count = T.get_node_attr(parent,
                                     'integer_infeasibility_count')
                integer_infeasibility_sum = T.get_node_attr(parent,
                                     'integer_infeasibility_sum')
                relax = T.get_node_attr(parent, 'lp_bound')
            elif status is 'fathomed':
                if T._incumbent_value is None:
                    print('WARNING: Encountered "fathom" line before'+\
                        ' first incumbent.')
                    print('  This may indicate an error in the input file.')
            elif status is 'integer':
                integer_infeasibility_count = None
                integer_infeasibility_sum = None
            T.AddOrUpdateNode(cur_index, parent, _direction[sense],
                                 status, relax,
                                 integer_infeasibility_count,
                                 integer_infeasibility_sum,
                                 branch_var = branch_var,
                                 branch_var_value = var_values[branch_var],
                                 sense = sense, rhs = rhs, obj = relax,
                                 color = color, style = 'filled',
                                 label = label, fillcolor = color)
            if status is 'integer':
                T._previous_incumbent_value = T._incumbent_value
                T._incumbent_value = relax
                T._incumbent_parent = parent
                T._new_integer_solution = True
            # Currently Broken
#           if ETREE_INSTALLED and T.attr['display'] == 'svg':
#               T.write_as_svg(filename = "node%d" % iter_count,
#                                 prevfile = "node%d" % (iter_count - 1),
#                                 nextfile = "node%d" % (iter_count + 1),
#                                 highlight = cur_index)
            if T.get_layout() == 'dot2tex':
                _dot2tex_label = {'>=':' \geq ', '<=':' \leq '}
                T.set_edge_attr(parent, cur_index, 'label',
                                   str(branch_var) + _dot2tex_label[sense] +
                                   str(rhs))
            else:
                T.set_edge_attr(parent, cur_index, 'label',
                                   str(branch_var) + sense + str(rhs))
        iter_count += 1
        if BBstatus == 'C':
            # Branching:
            # Choose a variable for branching
            branching_var = None
            if branch_strategy == FIXED_BRANCHING:
                #fixed order
                for i in VARIABLES:
                    frac = min(var[i].varValue-math.floor(var[i].varValue),
                               math.ceil(var[i].varValue) - var[i].varValue)
                    if (frac > 0):
                        min_frac = frac
                        branching_var = i
                        # TODO(aykut): understand this break
                        break
            elif branch_strategy == MOST_FRACTIONAL:
                #most fractional variable
                min_frac = -1
                for i in VARIABLES:
                    frac = min(var[i].varValue-math.floor(var[i].varValue),
                               math.ceil(var[i].varValue)- var[i].varValue)
                    if (frac> min_frac):
                        min_frac = frac
                        branching_var = i
            elif branch_strategy == PSEUDOCOST_BRANCHING:
                scores = {}
                for i in VARIABLES:
                    # find the fractional solutions
                    if (var[i].varValue - math.floor(var[i].varValue)) != 0:
                        scores[i] = min(pseudo_u[i][0]*(1-var[i].varValue),
                                        pseudo_d[i][0]*var[i].varValue)
                    # sort the dictionary by value
                branching_var = sorted(list(scores.items()),
                                       key=lambda x : x[1])[-1][0]
            else:
                print("Unknown branching strategy %s" %branch_strategy)
                exit()
            if branching_var is not None:
                print("Branching on variable %s" %branching_var)
            #Create new nodes
            if search_strategy == DEPTH_FIRST:
                priority = (-cur_depth - 1, -cur_depth - 1)
            elif search_strategy == BEST_FIRST:
                priority = (-relax, -relax)
            elif search_strategy == BEST_ESTIMATE:
                priority = (-relax - pseudo_d[branching_var][0]*\
                                 (math.floor(var[branching_var].varValue) -\
                                      var[branching_var].varValue),
                            -relax + pseudo_u[branching_var][0]*\
                                 (math.ceil(var[branching_var].varValue) -\
                                      var[branching_var].varValue))
            node_count += 1
            Q.push(node_count, priority[0], (node_count, cur_index, relax, branching_var,
                    var_values[branching_var],
                    '<=', math.floor(var[branching_var].varValue)))
            node_count += 1
            Q.push(node_count, priority[1], (node_count, cur_index, relax, branching_var,
                    var_values[branching_var],
                    '>=', math.ceil(var[branching_var].varValue)))
            T.set_node_attr(cur_index, color, 'green')
        if T.root is not None and display_interval is not None and\
                iter_count%display_interval == 0:
            T.display(count=iter_count)

    timer = int(math.ceil((time.time()-timer)*1000))
    print("")
    print("===========================================")
    print("Branch and bound completed in %sms" %timer)
    print("Strategy: %s" %branch_strategy)
    if complete_enumeration:
        print("Complete enumeration")
    print("%s nodes visited " %node_count)
    print("%s LP's solved" %lp_count)
    print("===========================================")
    print("Optimal solution")
    #print optimal solution
    for i in sorted(VARIABLES):
        if opt[i] > 0:
            print("%s = %s" %(i, opt[i]))
    print("Objective function value")
    print(LB)
    print("===========================================")
    if T.attr['display'] is not 'off':
        T.display(count=iter_count)
    T._lp_count = lp_count
    return opt, LB

if __name__ == '__main__':    
    T = BBTree()
    #T.set_layout('dot2tex')
    #T.set_display_mode('file')
    T.set_display_mode('xdot')
    #T.set_display_mode('pygame')
    CONSTRAINTS, VARIABLES, OBJ, MAT, RHS = GenerateRandomMIP(rand_seed = 120)
    BranchAndBound(T, CONSTRAINTS, VARIABLES, OBJ, MAT, RHS,
                   branch_strategy = MOST_FRACTIONAL,
                   search_strategy = BEST_FIRST,
                   display_interval = 10000)

