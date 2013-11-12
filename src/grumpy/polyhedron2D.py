import numpy as np
from polyhedron import Vrep, Hrep
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from math import ceil, floor

class Polyhedron2D:
    def __init__(self, points = None, rays = None, A = None, b = None):
        have_rep = False
        if points is not None and rays is not None:
            self.vrep = Vrep(points, rays)
            self.hrep = Hrep(self.vrep.A, self.vrep.b)
            have_rep = True
        if A is not None and b is not None:
            self.hrep = Hrep(A, b)
            have_rep = True
        if not have_rep:
            raise RuntimeError('Must provide either generators or constraint',
                               'matrix')
        self.ray_indices = []
        for i in range(len(self.hrep.generators)):
            if not self.hrep.is_vertex[i]:
                self.ray_indices.append(i)
        self.min_point = None
        self.max_point = None
        self.plot_max = None
        self.plot_min = None

    def make_integer_hull(self):
        if self.min_point is None or self.max_point is None:
            self.determine_hull_size()
        v = []
        for i in range(int(ceil(self.min_point[0])), 
                       int(ceil(self.max_point[0]))+1):
            for j in range(int(ceil(self.min_point[1])), 
                           int(ceil(self.max_point[1]))+1):
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

    def determine_plot_size(self, padding = 1):
        if self.min_point is None or self.max_point is None:
            self.determine_hull_size()
        self.plot_max = self.max_point[:]
        self.plot_min = self.min_point[:]
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
            self.plot_max = np.maximum(self.max_point, vertex + lim*ray)
            self.plot_min = np.minimum(self.min_point, vertex + lim*ray)
        else:
            current_point = 0
            next_point = self.hrep.adj[0][0]
        for i in range(len(self.hrep.generators)):
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
                self.plot_max = np.maximum(self.max_point, vertex + lim*ray)
                self.plot_min = np.minimum(self.min_point, vertex + lim*ray)
                break
            if self.hrep.adj[current_point][0] != prev_point:
                next_point = self.hrep.adj[current_point][0]
            else:
                next_point = self.hrep.adj[current_point][1]
    
        self.plot_min = np.array([floor(self.plot_min[i]) - padding 
                                  for i in [0, 1]])
        self.plot_max = np.array([ceil(self.plot_max[i]) + padding 
                                  for i in [0, 1]])

    def draw(self, ax, color = 'blue', linestyle = 'solid'):
        if self.plot_max == None or self.plot_min == None:
            self.determine_plot_size()
        x, y = [], []
        if len(self.ray_indices) > 0:
            current_point = self.ray_indices[0]
            ray = self.hrep.generators[current_point]
            if self.hrep.is_vertex[self.hrep.adj[self.ray_indices[0]][0]]:
                vertex = self.hrep.generators[
                    self.hrep.adj[self.ray_indices[0]][0]]
                next_point = self.hrep.adj[self.ray_indices[0]][0]
            else:
                vertex = self.hrep.generators[
                    self.hrep.adj[self.ray_indices[0]][1]]
                next_point = self.hrep.adj[self.ray_indices[0]][1]
            if ray[0] < 0:
                x_lim = (vertex[0] - self.plot_min[0])/-float(ray[0])
            elif ray[0] > 0:
                x_lim = (self.plot_max[0] - vertex[0])/float(ray[0])
            else:
                x_lim = 10000
            if ray[1] < 0:
                y_lim = (vertex[1] - self.plot_min[1])/-float(ray[1])
            elif ray[1] > 0:
                y_lim = (self.plot_max[1] - vertex[1])/float(ray[1]) 
            else:
                y_lim = 10000
            lim = min(x_lim, y_lim)*0.95
            ax.arrow(vertex[0], vertex[1], lim*ray[0], lim*ray[1], 
                     head_width = 0.02*(self.plot_max[0] - self.plot_min[0]), 
                     head_length = 0.03*(self.plot_max[0] - self.plot_min[0]), 
                     color = color, linestyle = linestyle)
        else:
            current_point = 0
            next_point = self.hrep.adj[0][0]
            x.append(self.hrep.generators[current_point][0])
            y.append(self.hrep.generators[current_point][1])
        for i in range(len(self.hrep.generators)):
            prev_point = current_point
            current_point = next_point
            if not self.hrep.is_vertex[current_point]:
                ray = self.hrep.generators[current_point]
                vertex = self.hrep.generators[prev_point]
                if ray[0] < 0:
                    x_lim = (vertex[0] - self.plot_min[0])/-float(ray[0])
                elif ray[0] > 0:
                    x_lim = (self.plot_max[0] - vertex[0])/float(ray[0])
                else:
                    x_lim = 10000
                if ray[1] < 0:
                    y_lim = (vertex[1] - self.plot_min[1])/-float(ray[1])
                elif ray[1] > 0:
                    y_lim = (self.plot_max[1] - vertex[1])/float(ray[1])
                else:
                    y_lim = 10000 
                lim = min(x_lim, y_lim)*0.95
                ax.arrow(vertex[0], vertex[1], lim*ray[0], lim*ray[1], 
                         head_width = 0.02*(self.plot_max[0] - 
                                            self.plot_min[0]), 
                         head_length = 0.03*(self.plot_max[0] - 
                                             self.plot_min[0]), 
                         color = color, linestyle = linestyle)
                break
            x.append(self.hrep.generators[current_point][0])
            y.append(self.hrep.generators[current_point][1])
            if self.hrep.adj[current_point][0] != prev_point:
                next_point = self.hrep.adj[current_point][0]
            else:
                next_point = self.hrep.adj[current_point][1]
        if linestyle == 'dashed':
            linestyle = '--'
        if linestyle == 'solid':
            linestyle = '-'
        line = lines.Line2D(x, y, color = color, linestyle = linestyle)
        ax.add_line(line)

def add_line(ax, coeffs, level, plot_max = None, plot_min = None, 
             color = 'blue', linestyle = 'solid'):
    if plot_max == None or plot_min == None:
        print 'Must have plot_max and plot_min set in order to add line'
        return
    x_intercept = [float(level - plot_min[1]*coeffs[1])/coeffs[0],
                   float(level - plot_max[1]*coeffs[1])/coeffs[0]]
    y_intercept = [float(level - plot_min[0]*coeffs[0])/coeffs[1],
                   float(level - plot_max[0]*coeffs[0])/coeffs[1]]
    
    x, y = [], []
    if coeffs[0]/coeffs[1] < 0:
        if y_intercept[1] > plot_max[1]:
            y.append(plot_max[1])
            x.append(float(level - plot_max[1]*coeffs[1])/coeffs[0])
        elif y_intercept[1] < plot_min[1]:
            return
        else:
            x.append(plot_max[0])
            y.append(float(level - plot_max[0]*coeffs[0])/coeffs[1])
        if y_intercept[0] < plot_min[1]:
            y.append(plot_min[1])
            x.append(float(level - plot_min[1]*coeffs[1])/coeffs[0])
        elif y_intercept[0] > plot_max[1]:
            return
        else:
            x.append(plot_min[0])
            y.append(float(level - plot_min[0]*coeffs[0])/coeffs[1])
    else:
        if y_intercept[1] < plot_min[1]:
            y.append(plot_min[1])
            x.append(float(level - plot_min[1]*coeffs[1])/coeffs[0])
        elif y_intercept[1] > plot_max[1]:
            return
        else:
            x.append(plot_max[0])
            y.append(float(level - plot_max[0]*coeffs[0])/coeffs[1])
        if y_intercept[1] > plot_max[1]:
            y.append(plot_max[1])
            x.append(float(level - plot_max[1]*coeffs[1])/coeffs[0])
        elif y_intercept[0] < plot_min[1]:
            return
        else:
            x.append(plot_min[0])
            y.append(float(level - plot_min[0]*coeffs[0])/coeffs[1])

    if linestyle == 'dashed':
        linestyle = '--'
    if linestyle == 'solid':
        linestyle = '-'
    line = lines.Line2D(x, y, color = color, linestyle = linestyle)
    ax.add_line(line)

if __name__ == '__main__':
    #points = np.random.random ((20,2))
    #p = make_hull(points)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    p = Polyhedron2D(A = [[4, 1], [1, 4], [1, -1], [-1, 0], [0, -1]], 
                   b = [28, 27, 1, 0, 0])
    p.draw(ax, color = 'blue', linestyle = 'solid')
    ax.set_xlim(p.plot_min[0], p.plot_max[0])
    ax.set_ylim(p.plot_min[1], p.plot_max[1])
    pI = p.make_integer_hull()
    pI.draw(ax, color = 'red', linestyle = 'dashed') 
    add_line(ax, [1, 1], 1, p.plot_max - [0.2, 0.2], p.plot_min + [0.2, 0.2], 
             linestyle = 'dashed')
    plt.show()


