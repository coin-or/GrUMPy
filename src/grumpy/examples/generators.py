import matplotlib.pyplot as plt
from src.grumpy.polyhedron2D import Polyhedron2D, add_line, add_point

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
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    p = Polyhedron2D(A = A, b = b)
    q = Polyhedron2D(points = points)
    p.draw(ax, color = 'blue', linestyle = 'solid', label = '$P$')
    q.draw(ax, color = 'black', linestyle = 'solid', color_int_points = True,
           label = '$Q$')
    ax.set_xlim(p.plot_min[0]-0.5, p.plot_max[0]+1)
    ax.set_ylim(p.plot_min[1], p.plot_max[1]+1)
    pI = p.make_integer_hull()
    pI.draw(ax, color = 'red', linestyle = 'dashed', label = 'conv(S)') 
    plt.text(-0.2, 1.3, r'$e_1 = (1.175, 1.375)$')
    plt.text(1.8, -0.2, r'$e_2 = (2, 0)$')
    plt.text(1.5, 2.5, r'$r_1 = (1, 1)$')
    plt.text(3.5, 0.5, r'$r_2 = (2, 1)$')
    add_point(ax, (1.175,1.375), 0.02, 'red')
    add_point(ax, (2,0), 0.02, 'red')
    plt.legend()
    plt.show()

