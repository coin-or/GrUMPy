import matplotlib.pyplot as plt
try:
    from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure
except ImportError:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure

if __name__ == '__main__':
    #points = np.random.random ((20,2))
    #p = make_hull(points)
    A = [[4, 1],
         [1, 4],
         [1, -1],
         [-1, 0],
         [0, -1]] 

    b = [28,
         27,
         1,
         1,
         1]
    f = Figure()
    p = Polyhedron2D(A = A, b = b)
    f.add_polyhedron(p, color = 'blue', linestyle = 'solid')
    f.set_xlim((2.5, 6.5))
    f.set_ylim((3, 7))
    #f.set_xlim(p.xlim)
    #f.set_ylim(p.ylim)
    pI = p.make_integer_hull()
    f.add_polyhedron(pI, color = 'red', linestyle = 'dashed',
                     show_int_points = True) 
    #f.add_point([5.666,5.333], 0.02, 'red')
    #f.add_text([5.7, 5.4], r'$(17/3, 16/3)$')

    #Lift and Project cut from infeasible basis 
    f.add_point((6.2,5.2), 0.02, 'red')
    f.add_line([4, 1], 28, linestyle = 'dashed')
    f.add_line_segment([5.6666, 5.333], [6.2, 5.2], linestyle = 'dashed')
    f.add_line_segment([5.7, 4.7], [6.2, 5.2], linestyle = 'dashed')
    f.add_line([1, 3], 21, p.xlim, p.ylim,
             color = 'black', linestyle = 'dashed')
    
    #This show the disjunction associated with row 2 to show how the cut in
    #that row is obtained as an intersection cut
    #f.add_line([1, -1], 0, p.plot_max, p.plot_min,
    #            color = 'black', linestyle = 'dashed')
    #f.add_line([1, -1], 1, p.plot_max, p.plot_min,
    #            color = 'black', linestyle = 'dashed')
    #f.add_text([4.8, 3.7], r'$x_1 - x_2 = 1$')
    #f.add_text([3.6, 4.2], r'$x_1 - x_2 = 0$')
    
    # CG cuts
    #f.add_line([4, 2], 33, p.plot_max, p.plot_min,
    #           color = 'green', linestyle = 'dashed')
    #f.add_text([5.4, 6], r'$4x_1 + 2x_2 \leq 33$')
    #f.add_line([3, 2], 27, p.plot_max, p.plot_min,
    #           color = 'green', linestyle = 'dashed')
    #f.add_text([5.1, 6], r'$3x_1+2x_2 \leq 27$')
    #f.add_line([1, 2], 16, p.plot_max, p.plot_min,
    #         color = 'green', linestyle = 'dashed')
    #f.add_text([4.5, 5.8], r'$x_1 + 2x_2 \leq 16$')
    
    # GMI cuts
    #f.add_line([12, 33], 234, p.plot_max, p.plot_min,
    #         color = 'black', linestyle = 'dashed')
    #f.add_text([5.9, 5], r'$12x_1 + 33x_2 \leq 234$')
    #f.add_line([3, 2], 27, p.plot_max, p.plot_min,
    #         color = 'black', linestyle = 'dashed')
    #f.add_line([3, 2], 26, p.plot_max, p.plot_min,
    #         color = 'black', linestyle = 'dashed')
    #f.add_text([4.7, 6], r'$3x_1 + 2x_2 \leq 26$')
    f.show()

    #x * 12 + y * 3 = 1
    #x * 33 + y * 2 = 2

    #y = 1 - x*33/2
    #12x + 3(1- 33x/2) = 1
    #2 = (3*33/2 - 12)x

    #x = 4/75
    #y = 9/75    

    
    

