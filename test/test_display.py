from grumpy import BBTree

if __name__=='__main__':
    bt = BBTree(display='pygame')
    CONSTRAINTS, VARIABLES, OBJ, MAT, RHS = bt.GenerateRandomMIP(rand_seed = 3)
    bt.BranchAndBound(CONSTRAINTS, VARIABLES, OBJ, MAT, RHS,
                     branch_strategy = 'Pseudocost', search_strategy = 'Best First',
                     display_interval = 0.1)

