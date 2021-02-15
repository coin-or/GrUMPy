import matplotlib.pyplot as plt
from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure

if __name__ == '__main__':
    #points = np.random.random ((20,2))
    #p = make_hull(points)
    A = [[-5, -3],
         [-5, 5],
         [1, -2]] 

    b = [-10,
         1,
         2]

    points = [[1.175, 1.375],
              [2,0],
              [2.175, 2.375],
              [3,1],
              [3.175, 2.375],
              [4,1]]
        
    f = Figure()
    p = Polyhedron2D(A = A, b = b)
    q = Polyhedron2D(points = points)
    
    f.add_polyhedron(p, color = 'blue', linestyle = 'solid', label = '$P$')
    f.add_polyhedron(q, color = 'black', linestyle = 'solid', show_int_points = True,
                     label = '$Q$')
    f.set_xlim([p.xlim[0], p.xlim[1]+1])
    f.set_ylim([p.ylim[0], p.ylim[1]+2])
    pI = p.make_integer_hull()
    f.add_polyhedron(pI, color = 'red', linestyle = 'dashed', label = 'conv(S)') 
    f.add_text((-0.2, 1.3), r'$e_1 = (1.175, 1.375)$')
    f.add_text((1.8, -0.2), r'$e_2 = (2, 0)$')
    f.add_text((1.5, 2.5), r'$r_1 = (1, 1)$')
    f.add_text((3.5, 0.5), r'$r_2 = (2, 1)$')
    f.add_point((1.175,1.375), 0.02, 'red')
    f.add_point((2,0), 0.02, 'red')
    f.show()

