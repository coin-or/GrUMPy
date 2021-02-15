from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure

from matplotlib import rcParams

from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

def make_legend_arrow(legend, orig_handle,
                      xdescent, ydescent,
                      width, height, fontsize):
    p = mpatches.FancyArrow(0, 0.5*height, width, 0,
                            length_includes_head=True,
                            head_width=0 )
    return p

rcParams["mathtext.fontset"] = 'cm'

A = [[ 1, -2],
     [-0, -1],
     [-2, -1],
     [-1, -0],
     [-1, -1],
     [1, -3],
     [0,  1],
     [1, 0]]

b = [4,
     0,   
     -3,
     0,
     -2,
     3,
     4,
     8]

c = [1, -0.6]
sense = ['Max', '<=']
integerIndices = [0, 1]

#This is for the separation code
#x0 = [1, 1]
x_sep = [0, 1]

points = [[1, 1],
                 [8, 2],
                 [3, 0],
                 [2, 0]
]

locations = [[0.5, 0.7],
             [8.1, 2.1],
             [2.8, -0.35],
             [1.6, -0.4]
]

d = [[0.2, -1.4],
     [-0.5, -1],
     [-1, -1]
     ]

P = Polyhedron2D(A = A, b = b)

iter = 3

f = Figure()

for iter in [1,2,3]:

    rays = [[0,0]]

    for i, e in enumerate(points):
        if (i <= iter):
            f.add_point(e, radius = 0.05, color = 'mediumblue')
            f.add_text(locations[i], f'$x^{i}$')
            if (i != 0):
                rays.append([e[0] - points[0][0],
                             e[1] - points[0][1]])
                
    C = f.ax.arrow(points[0][0], points[0][1], c[0], c[1],
                   head_width = 0.15, head_length = 0.15,
                   overhang = 0.3, color = 'navy', alpha = 0.7, linewidth = 2)
    
    D = f.ax.arrow(points[0][0], points[0][1],
                   d[iter-1][0], d[iter-1][1],
                   head_width = 0.15, head_length = 0.15,
                   overhang = 0.3, color = 'indianred',
                   alpha = 0.7, linewidth = 2)
    
    a = f.ax.arrow(points[0][0], points[0][1],
                   1.2*rays[1][0], 1.2*rays[1][1],
                   head_width = 0.15, head_length = 0.15,
                   overhang = 0.3, color = 'coral',
                   alpha = 0.7, linewidth = 2)
    #linestyle = (0, (1,2))
    
    if (iter == 2):
        f.ax.arrow(points[0][0], points[0][1],
                   1.2*rays[2][0], 1.2*rays[2][1],
                   head_width = 0.15, head_length = 0.15,
                   overhang = 0.3, color = 'coral',
                   alpha = 0.7, linewidth = 2)

    if (iter == 3):
        f.ax.arrow(points[0][0], points[0][1],
                   1.3*rays[3][0], 1.3*rays[3][1],
                   head_width = 0.15, head_length = 0.15,
                   overhang = 0.3, color = 'coral',
                   alpha = 0.7, linewidth = 2)
        
    #f.ax.annotate('', xytext = (points[0][0], points[0][1]),
    #              xy = (points[0][0]+1.4*rays[2][0],
    #                    points[0][1]+1.4*rays[2][1]), 
    #              arrowprops = {'arrowstyle':'->', 'linestyle':'--'},
    #               head_width = 0.2, head_length = 0.2, 
    #               color = 'red', linestyle = (5, (1,10))
    #)
    
    p = f.add_polyhedron(P, show_int_points = True,
                         color = 'steelblue', linewidth = 2)
    if (iter > 1):
        f.ax.fill([(points[i][0] + 0.4*rays[i][0]) for i in [0,1,iter]],
                  [(points[i][1] + 0.4*rays[i][1]) for i in [0,1,iter]],
                  color = 'peachpuff')
    
    f.set_xlim(P.xlim)
    f.set_ylim(P.ylim)

    #red_patch = mpatches.Patch(color='red', label='The red data', lw = 0.2,
    #                           linestyle = '-')
    #plt.legend(handles=[red_patch])
    
    if (iter == 1):
        plt.legend([p, a, C, D], ['$\mathrm{conv}(\mathcal{S}$)',
                                  '$\mathcal{C}_1(x^0)$', 
                                  '$c = d^0$', '$d^1$'],
                   handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow)})
    else:
        plt.legend([p, a, C, D], ['$\mathrm{conv}(\mathcal{S}$)',
                                  '$\mathcal{C}_'+f'{iter}(x^0)$', 
                                  '$c$', f'$d^{iter}$'],
                   handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow)})
        
    plt.savefig(f'inv-iter{iter}-new.png')
    f.show()
