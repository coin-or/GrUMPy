#!/bin/bash
python pprof.py -l 2 -c 4 --legend fixed_b_bestestimate\
                                   fixed_b_bestfirst\
                                   fixed_b_depthfirst\
                                   most_fractional_bestestimate\
                                   most_fractional_bestfirst\
                                   most_fractional_depthfirst\
                                   pseudocost_b_bestestimate\
                                   pseudocost_b_bestfirst\
                                   pseudocost_b_depthfirst\
                                   > num_lp_pp.eps
