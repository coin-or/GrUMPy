import numpy as np
from polyhedron import Vrep, Hrep
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from math import ceil, floor

class Polyhedron2D:
    def __init__(self, points = None, rays = None, A = None, b = None):
        have_rep = False
        if points is not None or rays is not None:
            self.vrep = Vrep(points, rays)
            self.hrep = Hrep(self.vrep.A, self.vrep.b)
            have_rep = True
        if A is not None and b is not None:
            self.hrep = Hrep(A, b)
            have_rep = True
        if not have_rep:
            raise RuntimeError('Must provide either generators or constraint'+
                               'matrix')
        self.ray_indices = []
        for i in range(len(self.hrep.generators)):
            if not self.hrep.is_vertex[i]:
                self.ray_indices.append(i)
        self.min_point = None
        self.max_point = None
        self.xlim = None
        self.ylim = None
        self.fig = None
        self.ax = None
        
    def make_integer_hull(self):
        if self.min_point is None or self.max_point is None:
            self.determine_hull_size()
        v = []
        for i in range(int(ceil(round(self.min_point[0], 3))), 
                       int(ceil(round(self.max_point[0], 3)))+1):
            for j in range(int(ceil(round(self.min_point[1], 3))), 
                           int(ceil(round(self.max_point[1], 3)))+1):
                if np.alltrue(np.dot(self.hrep.A, [i, j]) <= self.hrep.b):
                    v.append([i, j])
        r = [self.hrep.generators[self.ray_indices[i]].tolist() 
             for i in range(len(self.ray_indices))]
        return Polyhedron2D(points = v, rays = r)

    def determine_hull_size(self):
        self.min_point = np.array([10000, 10000])
        self.max_point = np.array([-10000, -10000])
        for i in range(len(self.hrep.generators)):
            if self.hrep.is_vertex[i]:
                self.max_point = np.maximum(self.max_point, self.hrep.generators[i])
                self.min_point = np.minimum(self.min_point, self.hrep.generators[i])

    def determine_plot_size(self, padding = [-1, 2]):
        if self.min_point is None or self.max_point is None:
            self.determine_hull_size()
        self.xlim = [self.min_point[0], self.max_point[0]]
        self.ylim = [self.min_point[1], self.max_point[1]]
        if len(self.ray_indices) > 0:
            current_point = self.ray_indices[0]
            ray = self.hrep.generators[current_point]
            if self.hrep.is_vertex[self.hrep.adj[self.ray_indices[0]][0]]:
                vertex = self.hrep.generators[self.hrep.adj[self.ray_indices[0]][0]]
                next_point = self.hrep.adj[self.ray_indices[0]][0]
            else:
                vertex = self.hrep.generators[
                    self.hrep.adj[self.ray_indices[0]][1]]
                next_point = self.hrep.adj[self.ray_indices[0]][1]
            if ray[0] < 0:
                x_lim = (vertex[0] - self.min_point[0])/-float(ray[0])
            elif ray[0] > 0:
                x_lim = (self.max_point[0] - vertex[0])/-float(ray[0])
            else:
                x_lim = 10000
            if ray[1] < 0:
                y_lim = (vertex[1] - self.min_point[1])/-float(ray[1])
            elif ray[1] > 0:
                y_lim = (self.max_point[1] - vertex[1])/-float(ray[1]) 
            else:
                y_lim = 10000
            lim = min(x_lim, y_lim) + 1
            plot_max = np.maximum(self.max_point, vertex + lim*ray)
            plot_min = np.minimum(self.min_point, vertex + lim*ray)
            self.xlim = [plot_min[0], plot_max[0]]
            self.ylim = [plot_min[1], plot_max[1]]
        else:
            current_point = 0
            next_point = self.hrep.adj[0][0]
        for _ in range(len(self.hrep.generators)):
            prev_point = current_point
            current_point = next_point
            if not self.hrep.is_vertex[current_point]:
                ray = self.hrep.generators[current_point]
                vertex = self.hrep.generators[prev_point]
                if ray[0] < 0:
                    x_lim = (vertex[0] - self.min_point[0])/-float(ray[0])
                elif ray[0] > 0:
                    x_lim = (self.max_point[0] - vertex[0])/float(ray[0])
                else:
                    x_lim = 10000
                if ray[1] < 0:
                    y_lim = (vertex[1] - self.min_point[1])/-float(ray[1])
                elif ray[1] > 0:
                    y_lim = (self.max_point[1] - vertex[1])/float(ray[1])
                else:
                    y_lim = 10000 
                lim = min(x_lim, y_lim) + 1
                plot_max = np.maximum(self.max_point, vertex + lim*ray)
                plot_min = np.minimum(self.min_point, vertex + lim*ray)
                self.xlim = [plot_min[0], plot_max[0]]
                self.ylim = [plot_min[1], plot_max[1]]
                break
            if self.hrep.adj[current_point][0] != prev_point:
                next_point = self.hrep.adj[current_point][0]
            else:
                next_point = self.hrep.adj[current_point][1]
    
        self.xlim = np.array(np.floor([self.xlim[i] + padding[i] for i in [0,1]]))
        self.ylim = np.array(np.floor([self.ylim[i] + padding[i] for i in [0,1]]))

class Figure:

    def __init__(self):
        self.fig = None

    def initialize(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.grid()

    def add_polyhedron(self, p, color = 'blue', linestyle = 'solid', label = None,
                       show_int_points = False):
        if self.fig is None:
            self.initialize()
        if p.xlim == None or p.ylim == None:
            p.determine_plot_size()
        x, y = [], []
        if len(p.ray_indices) > 0:
            current_point = p.ray_indices[0]
            ray = p.hrep.generators[current_point]
            if p.hrep.is_vertex[p.hrep.adj[p.ray_indices[0]][0]]:
                vertex = p.hrep.generators[
                    p.hrep.adj[p.ray_indices[0]][0]]
                next_point = p.hrep.adj[p.ray_indices[0]][0]
            else:
                vertex = p.hrep.generators[
                    p.hrep.adj[p.ray_indices[0]][1]]
                next_point = p.hrep.adj[p.ray_indices[0]][1]
            if ray[0] < 0:
                x_lim = (vertex[0] - p.xlim[0])/-float(ray[0])
            elif ray[0] > 0:
                x_lim = (p.xlim[1] - vertex[0])/float(ray[0])
            else:
                x_lim = 10000
            if ray[1] < 0:
                y_lim = (vertex[1] - p.ylim[0])/-float(ray[1])
            elif ray[1] > 0:
                y_lim = (p.ylim[1] - vertex[1])/float(ray[1]) 
            else:
                y_lim = 10000
            lim = min(x_lim, y_lim)*0.95
            self.ax.arrow(vertex[0], vertex[1], lim*ray[0], lim*ray[1], 
                          head_width = 0.02*(p.xlim[1] - p.xlim[0]), 
                          head_length = 0.03*(p.xlim[1] - p.xlim[0]), 
                          color = color, linestyle = linestyle)
        else:
            current_point = 0
            next_point = p.hrep.adj[0][0]
            x.append(p.hrep.generators[current_point][0])
            y.append(p.hrep.generators[current_point][1])
        for i in range(len(p.hrep.generators)):
            prev_point = current_point
            current_point = next_point
            if not p.hrep.is_vertex[current_point]:
                ray = p.hrep.generators[current_point]
                vertex = p.hrep.generators[prev_point]
                if ray[0] < 0:
                    x_lim = (vertex[0] - p.xlim[0])/-float(ray[0])
                elif ray[0] > 0:
                    x_lim = (p.xlim[1] - vertex[0])/float(ray[0])
                else:
                    x_lim = 10000
                if ray[1] < 0:
                    y_lim = (vertex[1] - p.ylim[0])/-float(ray[1])
                elif ray[1] > 0:
                    y_lim = (p.ylim[1] - vertex[1])/float(ray[1])
                else:
                    y_lim = 10000 
                lim = min(x_lim, y_lim)*0.95
                self.ax.arrow(vertex[0], vertex[1], lim*ray[0], lim*ray[1], 
                              head_width = 0.02*(p.ylim[1] - 
                                                 p.xlim[0]), 
                              head_length = 0.03*(p.ylim[1] - 
                                                  p.xlim[0]), 
                              color = color, linestyle = linestyle)
                break
            x.append(p.hrep.generators[current_point][0])
            y.append(p.hrep.generators[current_point][1])
            if p.hrep.adj[current_point][0] != prev_point:
                next_point = p.hrep.adj[current_point][0]
            else:
                next_point = p.hrep.adj[current_point][1]
        if linestyle == 'dashed':
            linestyle = '--'
        if linestyle == 'solid':
            linestyle = '-'
        line = lines.Line2D(x, y, color = color, linestyle = linestyle,
                            label = label)
        self.ax.add_line(line)
        if show_int_points:
            for i in range(int(ceil(p.min_point[0])), 
                           int(ceil(p.max_point[0]))+1):
                for j in range(int(ceil(p.min_point[1])), 
                               int(ceil(p.max_point[1]))+1):
                    if np.alltrue(np.dot(p.hrep.A, [i, j]) <= p.hrep.b):
                        self.add_point((i, j), radius = 0.02, color = 'black')

    def add_line_segment(self, point1, point2, color = 'blue', linestyle = 'solid',
                         label = None):
        if linestyle == 'dashed':
            linestyle = '--'
        if linestyle == 'solid':
            linestyle = '-'
        line = lines.Line2D([point1[0], point2[0]], 
                            [point1[1], point2[1]], 
                            color = color, linestyle = linestyle,
                            label = label)
        self.ax.add_line(line)
        
    def add_line(self, coeffs, level, xlim = None, ylim = None, 
                 color = 'blue', linestyle = 'solid', label = None):
        if self.fig is None:
            self.initialize()
        if xlim == None or ylim == None:
            print 'Must have plot_max and plot_min set in order to add line'
            return
        x_intercept = None
        y_intercept = None
        x = []
        y = []
        if coeffs[0] == 0 and coeffs[1] == 0:
            print 'Trying to plot line with zero coefficients...'
            return
        if coeffs[0] == 0:
            x = [xlim[0], xlim[1]]
            y = [level/coeffs[1], level/coeffs[1]]
        else:
            x_intercept = [float(level - ylim[0]*coeffs[1])/coeffs[0],
                           float(level - ylim[1]*coeffs[1])/coeffs[0]]
        if coeffs[1] == 0:
            x = [level/coeffs[0], level/coeffs[0]]
            y = [ylim[0], ylim[1]]
        else:
            y_intercept = [float(level - xlim[0]*coeffs[0])/coeffs[1],
                           float(level - xlim[1]*coeffs[0])/coeffs[1]]
        
        if x_intercept is not None and y_intercept is not None:
            if coeffs[0]/coeffs[1] < 0:
                if y_intercept[1] > ylim[1]:
                    y.append(ylim[1])
                    x.append(float(level - ylim[1]*coeffs[1])/coeffs[0])
                elif y_intercept[1] < ylim[0]:
                    return
                else:
                    x.append(xlim[1])
                    y.append(float(level - xlim[1]*coeffs[0])/coeffs[1])
                if y_intercept[0] < ylim[0]:
                    y.append(ylim[0])
                    x.append(float(level - ylim[0]*coeffs[1])/coeffs[0])
                elif y_intercept[0] > ylim[1]:
                    return
                else:
                    x.append(xlim[0])
                    y.append(float(level - xlim[0]*coeffs[0])/coeffs[1])
            else:
                if y_intercept[1] < ylim[0]:
                    y.append(ylim[0])
                    x.append(float(level - ylim[0]*coeffs[1])/coeffs[0])
                elif y_intercept[1] > ylim[1]:
                    return
                else:
                    x.append(xlim[1])
                    y.append(float(level - xlim[1]*coeffs[0])/coeffs[1])
                if y_intercept[1] > ylim[1]:
                    y.append(ylim[1])
                    x.append(float(level - ylim[1]*coeffs[1])/coeffs[0])
                elif y_intercept[0] < ylim[0]:
                    return
                else:
                    x.append(xlim[0])
                    y.append(float(level - xlim[0]*coeffs[0])/coeffs[1])
    
        if linestyle == 'dashed':
            linestyle = '--'
        if linestyle == 'solid':
            linestyle = '-'
        line = lines.Line2D(x, y, color = color, linestyle = linestyle,
                            label = label)
        self.ax.add_line(line)
    
    def add_point(self, center, radius = .02, color = 'red', ):
        if self.fig is None:
            self.initialize()
        self.ax.add_patch(plt.Circle(center, radius = radius, color = color))
    
    def add_text(self, loc, text):
        plt.text(loc[0], loc[1], text)
        
    def set_xlim(self, xlim):
        self.ax.set_xlim(xlim)

    def set_ylim(self, ylim):
        self.ax.set_ylim(ylim)

    def show(self):
        plt.legend()
        plt.show()
        self.fig = None

if __name__ == '__main__':
    #points = np.random.random ((20,2))
    #p = make_hull(points)

    f = Figure()
    p = Polyhedron2D(A = [[4, 1], [1, 4], [1, -1], [-1, 0], [0, -1]], 
                     b = [28, 27, 1, 0, 0])
    f.add_polyhedron(p, color = 'blue', linestyle = 'solid', label = 'p',
                     show_int_points = True)
    f.set_xlim(p.xlim)
    f.set_ylim(p.ylim)
    pI = p.make_integer_hull()
    f.add_polyhedron(pI, color = 'red', linestyle = 'dashed', label = 'pI') 
    f.add_point((5.666,5.333), 0.02, 'red')
    f.add_text((5.7, 5.4), r'$(17/3, 16/3)$')
    f.add_line([3, 2], 27, p.xlim, p.ylim,
               color = 'green', linestyle = 'dashed')
    f.show()
