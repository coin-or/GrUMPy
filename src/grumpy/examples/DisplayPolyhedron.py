'''
Created on Apr 4, 2019

@author: tkral
'''

from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure
points = [[2.5, 4.5], [6.5, 0.5], [0.5, 1],
          [7, 5.7], [7.7, 5], [2, 0.25]]
rays = []
c = [2, 5]
opt = [7, 5.7]
loc = (opt[0]+0.1, opt[1]-0.1)
obj_val = 42.5

p = Polyhedron2D(points = points, rays = rays)
f = Figure()
f.add_polyhedron(p, label = 'Polyhedron $P$', color = 'red')
f.set_xlim([p.xlim[0], p.xlim[1]+1])
f.set_ylim([p.ylim[0], p.ylim[1]+2])
f.add_line(c, obj_val, p.xlim + [0.2, 0.8], p.ylim + [0.2, 1.8], 
           linestyle = 'dashed', color = 'black', label = "Objective Function")
f.add_point(opt, 0.04, 'red')
f.add_text(loc, r'$x^* = %s$' % str(opt))
f.show()