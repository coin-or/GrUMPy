from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import hex
from builtins import range
from past.utils import old_div
# BBTree.py
#
#  Copyright 2009, 2010 Google Inc.
#  Copyright 2007 University of Pittsburgh.
#  Copyright 2012, 2013 Lehigh University
#  Google coding done by Brady Hunsaker.
#  U of Pittsburgh coding done by Osman Ozaltin and Brady Hunsaker.
#  Lehigh University coding done by Ted Ralphs and Aykut Bulut
#
#  This file was part of BAK (Branch-and-bound Analysis Kit).
#  It has now been incporporated into GrUMPy (Graphics for Understanding
#  Mathematical Programming in Python).
#
#  The contents of this file are subject to the Eclipse Public License
#  1.0.  (the "License"); you may not use this file except in
#  compliance with the License. You should have received a copy of
#  the Common Public License along with STOP.
#
#  Software distributed under the License is distributed on an "AS
#  IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or
#  implied. See the License for the specific language governing
#  rights and limitations under the License.
#
# For developers: Please keep the code style consistent with the Python
# style guide for Google's Summer of Code, except use 4 spaces to indent:
#   http://code.google.com/p/soc/wiki/PythonStyleGuide

__author__ = 'Brady Hunsaker, Osman Ozaltin, Ted Ralphs, Aykut Bulut'
__maintainer__ = 'Aykut Bulut (aykut@lehigh.edu)'

"""
This package is for visualizing branch-and-bound. It also contains
a basic branch-and-bound implementation primarily for classroom and educational
use.

Communication with solvers is through a grammar described in separate
documentation. Solvers can interface to this class in a number of
different ways and a number of different types of images may be created.

Images at intervals that can be specified on the command line as well as after
new incumbent solutions are found.

Note that the generation of tree images takes significantly longer than other
images because every node appears in the image.
"""

import math
import random
import re
import subprocess
import sys, os
import tempfile
import optparse
try:
    from src.gimpy import BinaryTree, XDOT_INSTALLED, MATPLOTLIB_INSTALLED, DOT2TEX_INSTALLED
    from src.gimpy import PIL_INSTALLED, ETREE_INSTALLED
    from src.gimpy import quote_if_necessary as quote
except ImportError:
    from coinor.gimpy import BinaryTree, XDOT_INSTALLED, MATPLOTLIB_INSTALLED, DOT2TEX_INSTALLED
    from coinor.gimpy import PIL_INSTALLED, ETREE_INSTALLED
    from coinor.gimpy import quote_if_necessary as quote
from io import StringIO
from pulp import LpVariable, lpSum, LpProblem, LpMaximize, LpConstraint
from pulp import LpStatus, value
from .forecasting import ForecastingChainedSequences

if MATPLOTLIB_INSTALLED:
    import matplotlib.pyplot as plt

if PIL_INSTALLED:
    from PIL import Image as PIL_Image

# branch strategy
BRANCH_STRATEGY = None
# search strategy
SEARCH_STRATEGY = None
# branching strategies
MOST_FRACTIONAL = 'Most Fraction'
FIXED_BRANCHING = 'Fixed Branching'
PSEUDOCOST_BRANCHING = 'Pseudocost Branching'
# search strategies
DEPTH_FIRST = 'Depth First'
BEST_FIRST = 'Best First'
BEST_ESTIMATE = 'Best Estimate'
INFINITY = sys.maxsize

class BBTree(BinaryTree):
    """
    Methods to process and visualize information about a b&b tree. It can
    process an output file (in a specific format, see BAK project in COIN-OR) of
    a solver that has three information. See run.py in examples directory fot
    this use. Moreover it implements a branch and bound method that can solve
    binary programs (0-1 variables only) using PuLP as an LP solver. It provides
    different branching and searching strategies. See test_strategies.py in test
    directory.

    This is the main class of GrUMPy. It inherits BinaryTree from GIMPy and
    keeps the entire branch-and-bound tree in self.
    """
    def __init__(self, **attrs):
        if not 'layout' in attrs:
            attrs['layout']= 'bak'
        BinaryTree.__init__(self, **attrs)
        # User-controlled constant values
        self._label = ''
        self._filename = None
        self._logscaley = False
        self._fathom = False
        self._wait_for_keypress = True
        self._edge_limit = 1000000
        # current time, updated each time we read a new line
        self._time = 0.0
        # use at most NUM nodes of each type in tree images; zero means no limit
        self._sample_tree = 0
        # Instance-dependent constant values
        self._optimization_sense = None
        self._start_time = None
        self._histogram_lower_bound = None
        self._histogram_upper_bound = None
        self._scatterplot_lower_bound = None
        self._scatterplot_upper_bound = None
        self._integer_infeasibility_lower_bound = None
        self._integer_infeasibility_upper_bound = None
        # Changing reference values
        self._image_counter = 0
        self._next_image_time = None
        self._incumbent_value = None
        self._incumbent_parent = None
        self._new_integer_solution = False
        self._max_objective_value = None
        self._min_objective_value = None
        self._max_integer_infeasibility_sum = None
        # List of incumbent path data files, for the all incumbent paths image
        self._incumbent_path_datafiles = []
        # Objects for measuring and predicting progress
        self._objective_gap_forecaster = ForecastingChainedSequences()
        self._sum_subtree_gaps_forecaster = ForecastingChainedSequences()
        self._sum_subtree_gaps_scale = 1.0
        self._previous_incumbent_value = None  # Only needed for SSG
        if 'display' in attrs:
            self.set_display_mode(attrs['display'])
        else:
            self.set_display_mode('off')

    def process_file(self, file_name):
        self._filename = file_name
        input_file = open(file_name, 'r')
        # Parse all the lines
        for line in input_file:
            self.ProcessLine(line)
            if self.root is not None:
                self.display()
        input_file.close()

    def write_as_dynamic_gexf(self, filename, mode = "Dot"):
        if not GEXF_INSTALLED:
            print('Gexf not installed. Exiting.')
            return
        if mode == 'Dot':
            try:
                gexf = Gexf("Mike O'Sullivan", "Dynamic graph file")
                graph = gexf.addGraph("directed", "dynamic", "Dynamic graph")
                objAtt = graph.addNodeAttribute("obj", "0.0", "float")
                currAtt = graph.addNodeAttribute("current", "1.0",
                                                 "integer", "dynamic")
                node_names = self.get_node_list()
                for name in node_names:
                    node = self.get_node(name)
                    if node.get("step") is None:
                        raise Exception("Node without step in BBTree",
                                        "node =", node)
                    curr_step = '%s' % node.get("step")
                    next_step = "%s" % (node.get("step") + 1)
                    n = graph.addNode(name, node.get_label(), start=curr_step)
                    if node.get("obj") is None:
                        raise Exception("Node without objective in BBTree",
                                        "node =", node)
                    n.addAttribute(objAtt, "%s" % node.get("obj"))
                    n.addAttribute(currAtt, "1", start=curr_step, end=next_step)
                    n.addAttribute(currAtt, "0", start=next_step)
                edge_names = self.get_edge_list()
                for i, (m_name, n_name) in enumerate(edge_names):
                    edge = self.get_edge(m_name, n_name)
                    if edge.get("step") is None:
                        raise Exception("Edge without step in BBTree",
                                        "edge =", (m_name, n_name))
                    curr_step = "%s" % edge.get("step")
                    graph.addEdge(i, m_name, n_name, start=curr_step)
                output_file = open(filename + ".gexf", "w")
                gexf.write(output_file)
            except Exception as e:
                print(e)
                print("No .gexf file created")
        else:
            raise Exception("Only Dot mode supported in write_as_dynamic_gexf")

    def set_display_mode(self, mode):
        if mode == 'off':
            self.attr['display'] = mode
        elif mode == 'matplotlib':
            if MATPLOTLIB_INSTALLED:
                self.attr['display'] = 'matplotlib'
            else:
                print('Matplotlib is not installed. Display is set to off.')
                self.attr['display'] = 'off'
        elif mode == 'PIL':
            if PIL_INSTALLED:
                self.attr['display'] = 'PIL'
            else:
                print('PIL is not installed. Display is set to off.')
                self.attr['display'] = 'off'
        elif mode == 'xdot':
            if XDOT_INSTALLED:
                self.attr['display'] = 'xdot'
            else:
                print('Xdot is not installed. Display is set to off.')
                self.attr['display'] = 'off'
        elif mode == 'file':
            self.attr['display'] = 'file'
        elif mode == 'matplotlib':
            self.attr['display'] = 'matplotlib'
        else:
            raise Exception('%s is not a valid display mode.' %mode)

    def display(self, item = 'all', basename = 'graph', format='png',
                count=None, pause=False, wait_for_click=True):
        '''
        Displays/Saves images requested. BranchAndBound method calls this method
        to visualize the branch and bound tree.
        '''
        if self.attr['layout'] != 'bak':
            if 'init_log_cond' in self.root.attr:
                max_log_cond = 0
                for n in list(self.nodes.values()):
                    if 'init_log_cond' in n.attr:
                        max_log_cond = max(n.attr['init_log_cond'], max_log_cond)
                for n in list(self.nodes.values()):
                    if 'init_log_cond' in n.attr:
                        log_begin = n.attr['init_log_cond']
                        log_end = n.attr['final_log_cond']
                        normalized_cond = (1-old_div(log_begin,max_log_cond))
                        color = str(hex(int(normalized_cond*256))[2:]) if normalized_cond >= .0625 else '0' + str(hex(int(normalized_cond*256))[2:])
                        n.attr['label'] = '%.0f \n %.0f' % (log_begin, log_end)
                        n.attr['color'] = '#' + color*3
                        n.attr['fillcolor'] = '#' + color*3
                        n.attr['style'] = 'filled'
                    else:
                        n.attr['label'] = ' '
            BinaryTree.display(self, pause = pause, wait_for_click = wait_for_click)
            return
        if self.attr['display'] == 'off':
            return
        if self.attr['display'] == 'matplotlib':
            gnuplot_script = None
            if item=='all':
                self.display_all()
            elif item=='tree':
                gnuplot_script = self.GenerateTreeImage()
            elif item=='scatterplot':
                gnuplot_script = self.GenerateScatterplot()
            elif item=='histogram':
                gnuplot_script = self.GenerateHistogram()
            elif item=='incumbent':
                gnuplot_script = self.GenerateIncumbentPath()
            elif item=='forecast':
                gnuplot_script = self.GenerateForecastImages()
            else:
                raise Exception('Unknown display() method argument %s' %item)
            if gnuplot_script is not None:
                self.display_image(gnuplot_script)
            # clean auxilary files.
            histogram_files = [f for f in os.listdir(".")
                               if f.startswith("histogram")]
            incumbent_files = [f for f in os.listdir(".")
                               if f.startswith("incumbentpath")]
            scatterplot_files = [f for f in os.listdir(".")
                                 if f.startswith("scatterplot")]
            t_fathomed_files = [f for f in os.listdir(".")
                                if f.startswith("tree_fathomed")]
            t_infeasible_files = [f for f in os.listdir(".")
                                if f.startswith("tree_infeasible")]
            t_pregnant_files = [f for f in os.listdir(".")
                                if f.startswith("tree_pregnant")]
            t_integer_files = [f for f in os.listdir(".")
                               if f.startswith("tree_integer")]
            t_branched_files = [f for f in os.listdir(".")
                               if f.startswith("tree_branched")]
            bak_filelist = (histogram_files + incumbent_files +
                            scatterplot_files + t_fathomed_files +
                            t_integer_files + t_branched_files +
                            t_infeasible_files + t_pregnant_files)
            for f in bak_filelist:
                os.remove(f)
        elif self.attr['display'] == 'xdot':
            if XDOT_INSTALLED:
                window = xdot.DotWindow()
                window.set_dotcode(self.to_string())
                window.connect('destroy', gtk.main_quit)
                gtk.main()
            else:
                print('Error: xdot not installed. Display disabled.')
                self.attr['display'] = 'off'
        elif self.attr['display'] == 'file':
            if count is not None:
                basename = basename + '_' + str(count)
            if self.attr['layout'] == 'dot2tex':
                if DOT2TEX_INSTALLED:
                    if format != 'pdf' or format != 'ps':
                        print("Dot2tex only supports pdf and ps formats,"+\
                            "falling back to pdf")
                        format = 'pdf'
                    self.set_layout('dot')
                    tex = dot2tex.dot2tex(self.to_string(), autosize=True,
                                          texmode = 'math',
                                          template = DOT2TEX_TEMPLATE)
                    f = open(basename+'.tex', 'w')
                    f.write(tex)
                    f.close()
                    subprocess.call(['latex', basename])
                    if format == 'ps':
                        subprocess.call(['dvips', basename])
                    elif format == 'pdf':
                        subprocess.call(['pdflatex', basename])
                    self.set_layout('dot2tex')
                    # clean auxilary files.
                    aux_filelist = [basename+'.tex', basename+'.log',
                                    basename+'.dvi', basename+'.aux']
                    for f in aux_filelist:
                        os.remove(f)
                else:
                    print("Dot2tex not installed, falling back to graphviz")
                    self.set_layout('dot')
                    self.write(basename+'.'+format, self.get_layout(), format)
            else:
                gnuplot_script = None
                if item=='all':
                    self.display_all()
                elif item=='tree':
                    gnuplot_script = self.GenerateTreeImage()
                elif item=='scatterplot':
                    gnuplot_script = self.GenerateScatterplot()
                elif item=='histogram':
                    gnuplot_script = self.GenerateHistogram()
                elif item=='incumbent':
                    gnuplot_script = self.GenerateIncumbentPath()
                elif item=='forecast':
                    gnuplot_script = self.GenerateForecastImages()
                else:
                    raise Exception('Unknown display() method argument %s' %item)
                if gnuplot_script is not None:
                    self.write_image(gnuplot_script, basename+'.'+format)
                # clean auxilary files.
                histogram_files = [f for f in os.listdir(".")
                                   if f.startswith("histogram")]
                incumbent_files = [f for f in os.listdir(".")
                                   if f.startswith("incumbentpath")]
                scatterplot_files = [f for f in os.listdir(".")
                                     if f.startswith("scatterplot")]
                t_fathomed_files = [f for f in os.listdir(".")
                                    if f.startswith("tree_fathomed")]
                t_infeasible_files = [f for f in os.listdir(".")
                                      if f.startswith("tree_infeasible")]
                t_pregnant_files = [f for f in os.listdir(".")
                                    if f.startswith("tree_pregnant")]
                t_integer_files = [f for f in os.listdir(".")
                                   if f.startswith("tree_integer")]
                t_branched_files = [f for f in os.listdir(".")
                                    if f.startswith("tree_branched")]
                bak_filelist = (histogram_files + incumbent_files +
                                scatterplot_files + t_fathomed_files +
                                t_integer_files + t_branched_files +
                                t_infeasible_files + t_pregnant_files)
                for f in bak_filelist:
                    os.remove(f)
        else:
            raise Exception('Unknown display mode %s' %self.attr['display'])

    def display_all(self):
        '''
        Assumes all the images have the same size.
        '''
        print ('This function is deprecated and no longer functions')
        return

        # Old source, just in case
        # if not (MATPLOTLIB_INSTALLED and PIL_INSTALLED:
        #     print('Matplotlib or Pillow not installed. Display disabled')
        #     return
        # tree = self.GenerateTreeImage()
        # scatterplot = self.GenerateScatterplot()
        # histogram = self.GenerateHistogram()
        # incumbent = self.GenerateIncumbentPath()
        # if tree is not None:
        #     imTree = StringIO(tree)
        #     pTree = pygame.image.load(imTree, 'png')
        #     sTree = pTree.get_size()
        #     rTree = pygame.Rect(0,0,sTree[0],sTree[1])
        # if scatterplot is not None:
        #     imScatterplot = StringIO(scatterplot)
        #     pScatterplot = pygame.image.load(imScatterplot, 'png')
        #     sScatterplot = pScatterplot.get_size()
        #     rScatterplot = pygame.Rect(sTree[0],0,sScatterplot[0],
        #                                sScatterplot[1])
        # if histogram is not None:
        #     imHistogram = StringIO(histogram)
        #     pHistogram = pygame.image.load(imHistogram, 'png')
        #     sHistogram = pHistogram.get_size()
        #     rHistogram = pygame.Rect(0,sTree[1],sHistogram[0],sHistogram[1])
        # if incumbent is not None:
        #     imIncumbent = StringIO(incumbent)
        #     pIncumbent = pygame.image.load(imIncumbent, 'png')
        #     sIncumbent = pIncumbent.get_size()
        #     rIncumbent = pygame.Rect(sTree[0],sTree[1],sIncumbent[0],
        #                              sIncumbent[1])
        # screen = pygame.display.set_mode((sTree[0]+sTree[0], sTree[1]+sTree[1]))
        # if tree is not None:
        #     screen.blit(pTree, rTree)
        # if scatterplot is not None:
        #     screen.blit(pScatterplot, rScatterplot)
        # if histogram is not None:
        #     screen.blit(pHistogram, rHistogram)
        # if incumbent is not None:
        #     screen.blit(pIncumbent, rIncumbent)
        # pygame.display.flip()
        # if self._wait_for_keypress:
        #     pause = True
        #     print("Press any key to continue to next image (ESCAPE to disable pausing)")
        # else:
        #     pause = False
        # while pause:
        #     e = pygame.event.poll()
        #     if e.type == pygame.KEYDOWN:
        #         keystate = pygame.key.get_pressed()
        #         if keystate[pygame.K_ESCAPE] != 0:
        #             self._wait_for_keypress = False
        #         break
        #     if e.type == pygame.QUIT:
        #         pause = False
        #         pygame.display.quit()
        #         # sys.exit() exits the whole program and I (aykut) guess it is
        #         # not appropriate here.
        #         #sys.exit()

    def write_image(self, gnuplot_script, filename):
        if not (MATPLOTLIB_INSTALLED and PIL_INSTALLED):
            print('Either matplotlib or Pillow is not installed. Display disabled')
            return
        try:
            p = subprocess.run(['gnuplot'], capture_output = True,
                               input = bytearray(gnuplot_script, 'utf8'))
        except OSError:
            print('''Gnuplot executable not found.
Gnuplot must be installed and in your search path.
After installation, ensure that the PATH variable is properly set.''')
            return
        p.check_returncode()

        if p.stderr:
            print(p.stderr)        

        if isinstance(filename, str):
            with open(filename, "w+b") as f:
                f.write(p.stdout)
        else:
            filename.write(p.stdout)

    def display_image(self, gnuplot_script, pause = False, wait_for_click = True):

        if not (PIL_INSTALLED and MATPLOTLIB_INSTALLED):
            print('Warning: Either matplotlib or Pillow is not installed. Cannot display.')
            return
        
        tmp_fd, tmp_name = tempfile.mkstemp()
        tmp_file = os.fdopen(tmp_fd, 'w+b')
        self.write_image(gnuplot_script, tmp_file)
        tmp_file.close()
        im = PIL_Image.open(tmp_name)
        plt.figure(1)
        plt.clf()
        plt.axis('off')
        plt.imshow(im, interpolation='bilinear' #resample=True
                   #extent = (0, 100, 0, 100)
        )
        if wait_for_click == True:
            plt.draw()
            try:
                if plt.waitforbuttonpress(timeout = 10000):
                    plt.close()
                    exit()
            except:
                exit()
        else:
            plt.show(block=pause)
        im.close()

    def set_label(self, label):
        self._label = label

    def set_logscaley(self, boolean):
        self._logscaley = boolean

    def set_fathom(self, boolean):
        self._fathom = boolean

    def set_edge_limit(self, limit):
        self._edge_limit = limit

    def set_sample_tree(self, number):
        self._sample_tree = number

    def AddOrUpdateNode(self, id, parent_id, branch_direction, status, lp_bound,
                        integer_infeasibility_count, integer_infeasibility_sum,
                        condition_begin = None, condition_end = None,
                        **attrs):
        '''
        This method designed to update nodes (in BAK) but we use it for
        updating/adding arcs. This is because of the tree data structure the
        authors adopted in BAK.
        We can divide these attributes such that some will belong to the edge
        parent_id->id and the others belong to the id node. The following shows
        whether the attribute belongs to edge or node.
        branch direction -> edge
        status -> node
        lp_bound -> node
        integer_infeasibility_count -> node
        integer_infeasibility_sum -> node
        parent_id -> node
        '''
        if (condition_begin is not None) and (condition_end is not None):
            #Figure out the color
            attrs['init_log_cond'] = math.log(condition_begin, 10)
            attrs['final_log_cond'] = math.log(condition_end, 10)
        if id in self.neighbors:
            # node already exists, update attributes
            self.set_node_attr(id, 'status', status)
            self.set_node_attr(id, 'lp_bound', lp_bound)
            self.set_node_attr(id, 'integer_infeasibility_count',
                               integer_infeasibility_count)
            self.set_node_attr(id, 'integer_infeasibility_sum',
                               integer_infeasibility_sum)
            if (condition_begin is not None) and (condition_end is not None):
                self.set_node_attr(id, 'init_log_cond',
                                   math.log(condition_begin, 10))
                self.set_node_attr(id, 'final_log_cond',
                                   math.log(condition_end, 10))
        elif self.root is None:
            self.add_root(id, status = status, lp_bound = lp_bound,
                          integer_infeasibility_count = integer_infeasibility_count,
                          integer_infeasibility_sum = integer_infeasibility_sum,
                          subtree_root = None, **attrs)
        elif parent_id is not None:
            if branch_direction == 'L':
                self.add_left_child(id, parent_id, status = status,
                    lp_bound = lp_bound,
                    integer_infeasibility_count = integer_infeasibility_count,
                    integer_infeasibility_sum = integer_infeasibility_sum,
                    subtree_root = None, **attrs)
            elif branch_direction == 'R':
                self.add_right_child(id, parent_id, status = status,
                    lp_bound = lp_bound,
                    integer_infeasibility_count = integer_infeasibility_count,
                    integer_infeasibility_sum = integer_infeasibility_sum,
                    subtree_root = None, **attrs)
        else:
            print('this should not happen.')
            raise Exception()
        if lp_bound is not None:
            self.UpdateObjectiveValueLimits(lp_bound)
            # Set optimization sense if not yet set
            if self._optimization_sense is None:
                if lp_bound < self.root.get_attr('lp_bound'):
                    self._optimization_sense = 'max'
                elif lp_bound > self.root.get_attr('lp_bound'):
                    self._optimization_sense = 'min'
            if self._optimization_sense == 'min' and lp_bound < self.root.get_attr('lp_bound'):
                print("Switching guess about objective sense to maximization based on bound change")
                self._optimization_sense = 'max'
            if self._optimization_sense == 'max' and lp_bound > self.root.get_attr('lp_bound'):
                print("Switching guess about objective sense to minimization based on bound change")
                self._optimization_sense = 'max'
        if integer_infeasibility_sum is not None:
            if (self._max_integer_infeasibility_sum is None or
                integer_infeasibility_sum >
                self._max_integer_infeasibility_sum):
                self._max_integer_infeasibility_sum = integer_infeasibility_sum

    def IsBetterThan(self, value1, value2):
        """
        Returns True if value1 is better than value2 as an objective value.
        This depends on the optimization sense of the instance.
        Args:
          value1: Float.
          value2: Float.
        Returns:
          True if value1 is better than value2 as an objective value.
        """
        if self._optimization_sense is None:
            print("Optimization sense is not set, assuming sense is miniminzation")
            self._optimization_sense = 'min'
        if self._optimization_sense == 'min':
            return value1 < value2
        else:
            return value1 > value2

    def IsBetterThanIncumbent(self, value):
        """
        Returns True if the passed value is better than current incumbent.
        Args:
          value: Float to use for comparison.
        Returns:
          True if the passed value is better than the current incumbent.
          'Better' is determined by the sense of optimization.
        """
        if self._incumbent_value is None:
            return True
        else:
            return self.IsBetterThan(value, self._incumbent_value)

    def UpdateObjectiveValueLimits(self, value):
        """Updates the min and max objective values if appropriate.
        Args:
          value: Float objective value.
        """
        if self._max_objective_value is None:
            self._max_objective_value = value
            self._min_objective_value = value
        else:
            if value > self._max_objective_value:
                self._max_objective_value = value
            if value < self._min_objective_value:
                self._min_objective_value = value

    def AddProgressMeasures(self):
        # No progress measures if there is no incumbent yet
        if self._incumbent_value is None:
            return
        # Store sum-of-subtree-gaps
        # We need to traverse all nodes unfortunately
        # TODO(bhunsaker): check whether we can just traverse active nodes
        active_node_count = 0
        subtree_bounds = {}
        new_integer_ssg = 0  # Only needed if this is a new integer solution
        for node_id in self.get_node_list():
            status = self.get_node_attr(node_id, 'status')
            if status == 'candidate' or status == 'pregnant':
                lp_bound = self.get_node_attr(node_id, 'lp_bound')
                subtree_root = self.get_node_attr(node_id, 'subtree_root')
                # Optional check for fathomed nodes.
                if (self._fathom and
                    not self.IsBetterThanIncumbent(lp_bound)):
                    continue
                active_node_count += 1
                if (subtree_root not in subtree_bounds or
                    self.IsBetterThan(lp_bound, subtree_bounds[subtree_root])):
                    subtree_bounds[subtree_root] = lp_bound
                if self._new_integer_solution:
                    self.set_node_attr(node_id, 'subtree_root', id)
                    new_integer_ssg += abs(self._incumbent_value - lp_bound)
        # If we have a new integer solution, we need to compute what
        # the measure would be with the previous integer solution for
        # scaling purposes.
        if (self._new_integer_solution and
            self._previous_incumbent_value is not None):
            reference_value = self._previous_incumbent_value
        else:
            reference_value = self._incumbent_value
        sum_subtree_gaps = 0
        for lp_bound in list(subtree_bounds.values()):
            sum_subtree_gaps += abs(reference_value - lp_bound)
        # Start a new sequence if a new integer solution was just found
        if self._new_integer_solution:
            if new_integer_ssg >= 1e-6:
                scale_factor = (old_div(float(sum_subtree_gaps),
                                float(new_integer_ssg)))
            else:
                scale_factor = 1.0
            self._sum_subtree_gaps_forecaster.StartNewSequence(scale_factor)
            # sum_subtree_gaps was based on the previous integer solution;
            # update it now
            sum_subtree_gaps = new_integer_ssg
        self._sum_subtree_gaps_forecaster.AddMeasure(self._time,
                                                     sum_subtree_gaps,
                                                     active_node_count,
                                                     len(self.get_node_list()))
        # Add objective gap measure.  Note that this relies on the
        # active_node_count computed above.
        if self._new_integer_solution:
            self._objective_gap_forecaster.StartNewSequence(1.0)
        if self._optimization_sense == 'min':
            obj_gap = self._incumbent_value - self._min_objective_value
        else:
            obj_gap = self._max_objective_value - self._incumbent_value
        self._objective_gap_forecaster.AddMeasure(self._time, obj_gap,
                                                  active_node_count,
                                                  len(self.get_node_list()))

    def GetImageCounterString(self):
        """
        Returns a string with the image counter.
        """
        return '%03d' % self._image_counter

    def WriteHistogramScript(self, num_bins, bin_width, max_bin_count,
                                 lp_bound, data_filename, output_file):
        """
        Write a Gnuplot script file to generate a histogram image.
        Args:
          num_bins: Integer number of bins for the histogram.
          bin_width: Float width of the bins in terms of objective values.
          max_bin_count: Integer number of the highest bin count.
          lp_bound: Float value of the current LP bound.
          data_filename: String name of the file; used for display purposes.
        """
        # TODO(bhunsaker): add checks for bin_width zero
        if self._incumbent_value is not None:
            incumbent_bin = 1 + ((self._incumbent_value -
                                  self._histogram_lower_bound) // bin_width)
            incumbent_x_coord = 0.5 + (old_div((self._incumbent_value -
                                        self._histogram_lower_bound),
                                       bin_width))
        lp_bound_bin = 1 + ((lp_bound - self._histogram_lower_bound) //
                            bin_width)
        lp_bound_x_coord = 0.5 + (old_div((lp_bound - self._histogram_lower_bound),
                                  bin_width))
        # TODO(bhunsaker): Ask Osman about adjust_xcoord option, which appears
        #    to put the vertical lines at the edge of bins rather than the
        #    true location.
        # Output the Gnuplot script to a file.
        script = ""
        # Set terminal for the output files.
        script += 'set terminal png notransparent size 480,360\n\n'
        # Make settings for the scatter plot.
        index_string = self.GetImageCounterString()
        output_filename = "histogram."+index_string+".png"
        if output_file:
            script += 'set output "%s"\n' % output_filename
        if self._filename is None:
            script += 'set title "Histogram of LP Bounds"\n'
        else:
            script += ('set title "Histogram of LP Bounds: %s, %s, %.2fs"\n'
                       % (self._filename, self._label, self._time))
        script += 'set xlabel "obj. value"\n'
        script += 'set ylabel "number of nodes"\n'
        if self._logscaley:
            script += 'set logscale y\n'
        script += 'set nokey\n'
        script += 'set tics scale 0.001\n'

        script += 'set xrange [0:%d+1]\n' % num_bins
        if self._logscaley:
            script += 'set yrange [1:%d*1.2]\n' % max_bin_count
        else:
            script += 'set yrange [0:%d*1.2]\n' % max_bin_count
        script += 'set xtics rotate by 90\n'
        # Mark tics for each bin.
        script += 'set xtics ('
        # TODO(bhunsaker): Consider putting this in a loop.
        x_values = ['"%0.2f" %0.2f' %
                    (self._histogram_lower_bound + i * bin_width, i + 0.5)
                    for i in range(num_bins + 1)]
        script += ', '.join(x_values) + ')\n'
        # Plot LP bound and incumbent tics.
        script += 'set x2tics ('
        script += '"%0.2f" %d' % (lp_bound, lp_bound_bin)
        if self._incumbent_value is not None:
            script += ', "%0.2f"%d)\n' % (self._incumbent_value,
                                                  incumbent_bin)
        else:
            script += ')\n'
        plot_parts = []
        # Plot the data points.
        plot_parts.append('\'%s\' with boxes fill solid 0.2' % data_filename)
        # Draw the vertical lp_bound and incumbent lines.
        script += 'set parametric\n'
        script += 'set trange [0:%d*1.5]\n' % max_bin_count
        plot_parts.append('%0.2f,t linetype 2' % lp_bound_x_coord)
        if self._incumbent_value is not None:
            plot_parts.append('%0.2f,t linetype 5' % incumbent_x_coord)
        script += 'plot %s\n' % ', '.join(plot_parts)
        script += 'unset parametric\n'
        script += 'show output\n'
        return script

    def AdjustHistogramEndBins(self, objective_list, num_bins, bin_width,
                               bin_counts, bin_centers, bin_widths):
        """
        Adjusts the two end bins if necessary to make them narrower.
        The two end bins may need to be narrower than the other bins so that
        they do not go past the current incumbent value on one end and the
        current lp bound on the other.  So that the histogram is still correct
        in areas, the height of these bins needs to be adjusted so that the
        area does not change.

        Note that there is likely to be some bias toward taller bins on
        the ends since they always have a point at one end of their width.  It
        may be more accurate visually to ignore or discount that one point when
        determining the bin height, but that is not currently done.

        Args:
          objective_list: List of float objective values.
          num_bins: Integer number of bins.
          bin_width: Float standard width of bins in terms of objective values.
          bin_counts: List of integer counts for each bin.
          bin_centers: List of float coordinates for the center of each bin.
          bin_widths: List of float widths for bins, allowing for individualized
            widths.
        """
        if self._optimization_sense == 'min':
            lp_bound = min(objective_list)
            lower_bound = lp_bound
            if self._incumbent_value is not None:
                upper_bound = self._incumbent_value
            else:
                upper_bound = self._histogram_upper_bound
        else:
            lp_bound = max(objective_list)
            upper_bound = lp_bound
            if self._incumbent_value is not None:
                lower_bound = self._incumbent_value
            else:
                lower_bound = self._histogram_lower_bound
        # The end bins may have unusual centers and widths
        highest_nonempty_bin = int((upper_bound -
                                    self._histogram_lower_bound) // bin_width)
        if (highest_nonempty_bin < num_bins and
            bin_counts[highest_nonempty_bin] > 0):
            highest_x_coord = 0.5 + (old_div((upper_bound -
                                      self._histogram_lower_bound), bin_width))
            highest_nonempty_bin_width, unused_int = math.modf(0.5 +
                                                               highest_x_coord)
            if highest_nonempty_bin_width == 0.0:
                highest_nonempty_bin_width = 1.0
            bin_widths[highest_nonempty_bin] = highest_nonempty_bin_width
            bin_centers[highest_nonempty_bin] = highest_x_coord - (
                old_div(highest_nonempty_bin_width, 2))
            # Scale the height appropriately
            bin_counts[highest_nonempty_bin] /= bin_widths[highest_nonempty_bin]
        lowest_nonempty_bin = int((lower_bound -
                                   self._histogram_lower_bound) // bin_width)
        if bin_counts[lowest_nonempty_bin] > 0:
            lowest_x_coord = 0.5 + (old_div((lower_bound -
                                     self._histogram_lower_bound), bin_width))
            lowest_nonempty_bin_excess, unused_int = math.modf(0.5 +
                                                               lowest_x_coord)
            bin_widths[lowest_nonempty_bin] = 1.0 - lowest_nonempty_bin_excess
            bin_centers[lowest_nonempty_bin] = (
                lowest_x_coord + old_div(bin_widths[lowest_nonempty_bin], 2))
            # Scale the height appropriately
            bin_counts[lowest_nonempty_bin] /= bin_widths[lowest_nonempty_bin]

    def GenerateHistogram(self, output_file = False):
        """
        Generate files necessary for a histogram image.
        Two files are necessary: a data file and a Gnuplot script file (which
        references the data file).
        Args:
          time: Float number of seconds since the start of optimization.
        """
        num_bins = 20
        # Compute the bin width and counts.
        objective_list = []
        for n in self.get_node_list():
            if (self.get_node_attr(n,'status') == 'candidate' or
                self.get_node_attr(n,'status') == 'pregnant'):
                lp_bound = self.get_node_attr(n,'lp_bound')
                if not self.IsBetterThanIncumbent(lp_bound):
                    continue
                objective_list.append(lp_bound)
        # TODO(aykut) added the following check, we need it since we generate
        # histograms real time
        # we can not generate histogram if we do not have upperl and lower
        #bounds
        if len(objective_list)==0 or self._incumbent_value is None:
            return None
        # The first time we create a histogram, set bounds for objective
        # values.
        # TODO(bhunsaker): Consider bounds; talk to Osman.
        if self._histogram_lower_bound is None:
            if self._optimization_sense == 'min':
                self._histogram_lower_bound = min(objective_list)
                if self._incumbent_value is not None:
                    self._histogram_upper_bound = self._incumbent_value
                else:
                    self._histogram_upper_bound = max(objective_list)
            else:
                self._histogram_upper_bound = max(objective_list)
                if self._incumbent_value is not None:
                    self._histogram_lower_bound = self._incumbent_value
                else:
                    self._histogram_lower_bound = min(objective_list)
        bin_width = old_div((self._histogram_upper_bound -
                     self._histogram_lower_bound), float(num_bins))
        bin_counts = [0.0 for i in range(num_bins)]
        for value in objective_list:
            bin = int(math.floor(old_div((value - self._histogram_lower_bound),
                                 bin_width)))
            # Special case for the largest value.
            if (value >= self._histogram_upper_bound and
                value < self._histogram_upper_bound + 1e-6):
                bin = num_bins - 1
            if bin < 0:
                return
            assert bin < num_bins, '%d (%f) !< %d (%f)' % (
                bin, value, num_bins, self._histogram_upper_bound)
            bin_counts[bin] += 1
        max_bin_count = max(bin_counts)
        bin_centers = [i + 1.0 for i in range(len(bin_counts))]
        bin_widths = [1.0 for i in range(len(bin_counts))]
        self.AdjustHistogramEndBins(objective_list, num_bins, bin_width,
                                    bin_counts, bin_centers, bin_widths)
        if self._optimization_sense == 'min':
            lp_bound = min(objective_list)
        else:
            lp_bound = max(objective_list)
        # Output the bin data to a file.
        index_string = self.GetImageCounterString()
        data_filename = 'histogram%s.dat' % index_string
        data_file = open(data_filename, 'w')
        for index in range(len(bin_counts)):
            data_file.write('%f %f %f\n' % (bin_centers[index],
                                            bin_counts[index],
                                            bin_widths[index]))
        data_file.close()
        histogram_script = self.WriteHistogramScript(num_bins, bin_width,
                           max_bin_count, lp_bound, data_filename, output_file)
        # TODO(bhunsaker): Temporary hack
        #   This allows the bounds to be reset until an incumbent is found.
        if self._incumbent_value is None:
            self._histogram_lower_bound = None
            self._histogram_upper_bound = None
        return histogram_script

    def GetImageObjectiveBounds(self, min_value, max_value):
        """
        Return min and max bounds to be used for images.
        Images should use bounds that are slightly wider than observed
        objective values.  Also, the special case of a single value must be
        handled.
        Args:
          min_value: Float minimum objective value.
          max_value: Float maximum objective value.
        Returns:
          A tuple of two float values (lower_bound, upper_bound).
        """
        obj_range = max_value - min_value
        if obj_range > 0:
            image_max_obj = max_value + 0.01 * obj_range
            image_min_obj = min_value - 0.01 * obj_range
        else:
            if max_value >= 0:
                image_max_obj = 1.01 * max_value
            else:
                image_max_obj = 0.99 * max_value
            if min_value >= 0:
                image_min_obj = 0.99 * min_value
            else:
                image_min_obj = 1.01 * min_value
        return (image_min_obj, image_max_obj)

    def WriteScatterplotScript(self, data_filename, output_file):
        """
        Write a Gnuplot script file to generate a scatterplot image.
        Args:
          data_filename: String name of the file; used for display purposes.
        """
        image_min_obj, image_max_obj = self.GetImageObjectiveBounds(
            self._scatterplot_lower_bound, self._scatterplot_upper_bound)
        index_string = self.GetImageCounterString()
        output_filename = "scatterplot."+index_string+".png"
        script = ""
        # Set terminal for the output files.
        script += 'set terminal png notransparent size 480,360\n\n'
        # Make settings for the scatter plot.
        if output_file:
            script += 'set output "%s"\n' % output_filename
        if self._filename is None:
            script += 'set title "Scatterplot"\n'
        else:
            script += ('set title "Scatterplot: %s, %s, %ds"\n' % (
                    self._filename, self._label, int(self._time)))
        script += 'set pointsize 0.8\n'
        script += 'set nokey\n'
        script += 'set xlabel \"sum of int. infeas.\"\n'
        script += 'set ylabel \"obj. value\"\n'
        script += ('set xrange [0:%0.6f+2]\n' %
                          self._max_integer_infeasibility_sum)
        script += ('set yrange [%0.6f:%0.6f]\n' % (image_min_obj,
                                                          image_max_obj))
        plot_parts = []
        # Plot the data points.
        plot_parts.append('\'%s\' with points pointtype 2 linetype 1' %
                          data_filename)
        # Also plot the incumbent line.
        if self._incumbent_value is not None:
            plot_parts.append('%0.6f linetype 2 linewidth 0.5' %
                              self._incumbent_value)
        # Plot the incumbent's parent if it's available.
        if self._incumbent_parent is not None:
            #incumbent_parent = self.get_node(self._incumbent_parent)
            plot_parts.append('"< echo %0.6f %0.6f" '
                              'with points pointtype 9 pointsize 1.2' %
                              (self.get_node_attr(self._incumbent_parent,
                                                  'integer_infeasibility_sum'),
                               self.get_node_attr(self._incumbent_parent,
                                                  'lp_bound')))
        script += 'plot %s\n' % ', '.join(plot_parts)
        script += 'show output\n'
        return script

    def GenerateScatterplot(self, output_file = False):
        """
        Generate files necessary for a scatterplot image.
        Two files are necessary: a data file and a Gnuplot script file (which
        references the data file).
        Args:
            output_file: if not given the gnuplot image will not be written
        to disk but returned (to be displayed in matplotlib window)
        """
        # Output data points.
        index_string = self.GetImageCounterString()
        data_filename = 'scatterplot%s.dat' % index_string
        data_file = open(data_filename, 'w')
        if self._scatterplot_lower_bound is None:
            bounds = []
        # Write objective values and integer infeasibility sum information
        # for candidate and pregnant nodes.
        for node in self.get_node_list():
            status = self.get_node_attr(node, 'status')
            lp_bound = self.get_node_attr(node,'lp_bound')
            if status == 'candidate' or status == 'pregnant':
                # Optional check for fathomed nodes.
                if (self._fathom and
                    not self.IsBetterThanIncumbent(lp_bound)):
                    continue
                data_file.write('%0.6f %0.6f\n' % (
                        self.get_node_attr(node, 'integer_infeasibility_sum'),
                        lp_bound))
                # Set the image objective bounds the first image.
                if self._scatterplot_lower_bound is None:
                    bounds.append(lp_bound)
        data_file.close()
        if self._scatterplot_lower_bound is None:
            if len(bounds) <= 1:
                return None
            self._scatterplot_lower_bound = min(bounds)
            self._scatterplot_upper_bound = max(bounds)
            # The incumbent overrides a bound if present.
            if self._incumbent_value is not None:
                if self._optimization_sense == 'min':
                    self._scatterplot_upper_bound = self._incumbent_value
                else:
                    self._scatterplot_lower_bound = self._incumbent_value
        scatterplot_script = self.WriteScatterplotScript(data_filename,
                                                         output_file)
        return scatterplot_script

    def WriteIncumbentPathScript(self, data_filename):
        """
        Write a Gnuplot script file to generate an incumbent path image.
        Args:
          data_filename: String name of the file; used for display purposes.
        """
        image_min_obj, image_max_obj = self.GetImageObjectiveBounds(
            self._scatterplot_lower_bound, self._scatterplot_upper_bound)
        script = ''
        # Set terminal for the output files.
        script += 'set terminal png notransparent size 480,360\n\n'
        if self._filename is None:
            script += 'set title "Incumbent path"\n'
        else:
            script += ('set title "Incumbent path (%s %.2fs %s)"\n' % (
                    self._filename, self._time, self._label))
        script += 'set pointsize 0.8\n'
        script += 'set nokey\n'
        script += 'set xlabel \"sum of int. infeas.\"\n'
        script += 'set ylabel \"obj. value\"\n'
        script += ('set xrange [0:%0.6f+2]\n' %
                          self._max_integer_infeasibility_sum)
        script += ('set yrange [%0.6f:%0.6f]\n' % (image_min_obj,
                                                          image_max_obj))
        # Plot the data points and connecting lines.
        script += ('plot \'%s\' with points pointtype 2, '
                          '\'%s\' with lines linetype 2\n' %
                          (data_filename, data_filename))
        script += 'show output\n'
        return script

    def WriteAllIncumbentPathsScript(self):
        """
        Return a Gnuplot script string to generate an incumbent path image.
        Args:
          data_filenames: List of string names of files.
        """
        data_filenames = self._incumbent_path_datafiles
        image_min_obj, image_max_obj = self.GetImageObjectiveBounds(
            self._scatterplot_lower_bound, self._scatterplot_upper_bound)
        script = ''
        # Set terminal for the output files.
        script += 'set terminal png notransparent size 480,360\n\n'

        # Make settings for the scatter plot.
        if self._filename is None:
            script += 'set title "Incumbent paths"\n'
        else:
            script += ('set title "Incumbent paths (%s %.2fs %s)"\n' % (
                    self._filename, self._time, self._label))
        script += 'set pointsize 0.8\n'
        script += 'set nokey\n'
        script += 'set xlabel \"sum of int. infeas.\"\n'
        script += 'set ylabel \"obj. value\"\n'
        script += ('set xrange [0:%0.6f+2]\n' %
                          self._max_integer_infeasibility_sum)
        script += ('set yrange [%0.6f:%0.6f]\n' % (image_min_obj,
                                                          image_max_obj))
        # Plot the data points and connecting lines.
        command_list = []
        for filename in data_filenames:
            command_list.append('\'%s\' with points pointtype 2, '
                                '\'%s\' with lines linetype 2' %
                                (filename, filename))
        script += 'plot %s\n' % ','.join(command_list)
        script += 'show output\n'
        return script

    def GenerateIncumbentPath(self):
        """
        Generate files necessary for an incumbent scatterplot path image.
        Two files are necessary: a data file and a Gnuplot script file (which
        references the data file).
        """
        if self._incumbent_parent is None:
            return
        if self._scatterplot_lower_bound is None:
            return
        if self._scatterplot_upper_bound is None:
            return
        index_string = self.GetImageCounterString()
        # Output data points.
        data_filename = 'incumbentpath%s.dat' % index_string
        data_file = open(data_filename, 'w')
        # Write objective values and integer infeasibility sum information
        # for ancestor nodes.
        data_file.write('0 %0.6f\n' % self._incumbent_value)
        parent = self._incumbent_parent
        # TODO(bhunsaker): I think the following assumes a unique value for the
        #   parent of the root.
        while parent != None:
            data_file.write('%0.6f %0.6f\n'
                            % (self.get_node_attr(parent,
                               'integer_infeasibility_sum'),
                             self.get_node_attr(parent, 'lp_bound')))
            parent = self.get_node_attr(parent, 'parent')
        data_file.close()
        self._incumbent_path_datafiles.append(data_filename)
        # Output the Gnuplot script to a file.
        path_script = self.WriteIncumbentPathScript(data_filename)
        return path_script

    def GenerateAllIncumbentPaths(self):
        """
        Generate file for a path image with all incumbent paths.
        Data files were previously generated for each incumbent.  This re-uses
        those files.
        """
        all_path_script = self.WriteAllIncumbentPathsScript()

    def WriteTreeScript(self, additional_lines = None):
        """
        Write a Gnuplot script file to generate a tree image.
        Args:
          additional_lines: String with additional lines to be added to the
            script file.
        """
        image_min_obj, image_max_obj = self.GetImageObjectiveBounds(
            self._min_objective_value, self._max_objective_value)
        data = ''
        data += 'set terminal png notransparent size 480,360\n'
        data += 'set output "%s"\n' % output_file
        data += 'set nokey\n'
        data += 'set autoscale\n'
        data += 'set tics scale 0.001\n'
        data += 'set pointsize 0.5\n'
        data += 'set xrange [-0.1:1.1]\n'
        data += 'set yrange [%0.6f:%0.6f]\n' % (image_max_obj,
                                                          image_min_obj)
        data += 'set format x ""\n'
        data += 'set ylabel "obj. value"\n'
        if self._filename is None:
            data += 'set title "B&B tree"\n'
        else:
            data += 'set title "B&B tree (%s %.2fs %s)"\n\n' % (
                self._filename, self._time, self._label)
        for line in additional_lines:
            data += line
        return data

    def GetTreeFixedHorizontalPositions(self):
        """
        Returns horizontal positions for all nodes based on fixed positions.
        Returns:
          Dictionary of float horizontal positions, keyed by node id.
        """
        # Statistics needed for horizontal positions.
        horizontal_lower_bound = dict.fromkeys(self.get_node_list(), 0.0)
        horizontal_upper_bound = dict.fromkeys(self.get_node_list(), 1.0)
        horizontal_positions = dict.fromkeys(self.get_node_list())
        horizontal_positions[self.root.name] = 0.5
        # sort node list
        node_id_list = sorted(self.get_node_list())
        node_id_list_int = list(int(n) for n in node_id_list)
        node_id_list_int = sorted(node_id_list_int)
        node_id_list = list(str(n) for n in node_id_list_int)
        for node_id in node_id_list:
            if node_id == self.root.name:
                continue
            parent_id = self.get_node_attr(node_id, 'parent')
            branch_direction = self.get_node_attr(node_id, 'direction')
            if branch_direction == 'R':
                horizontal_lower_bound[node_id] = horizontal_positions[
                    parent_id]
                horizontal_upper_bound[node_id] = horizontal_upper_bound[
                    parent_id]
            elif branch_direction == 'L':
                horizontal_lower_bound[node_id] = horizontal_lower_bound[
                    parent_id]
                horizontal_upper_bound[node_id] = horizontal_positions[
                    parent_id]
            else:
                print('Error: node %s has unsupported branching direction.'\
                    %node_id)
                print('Fixed-position tree images only support L and R ')
                print('branching.')
                sys.exit(1)
            horizontal_positions[node_id] = old_div((
                horizontal_upper_bound[node_id] +
                horizontal_lower_bound[node_id]), 2)
        return horizontal_positions

    def GetTreeHorizontalPositions(self):
        """
        Returns horizontal positions for all nodes.
        Each node is given equal horizontal space.
        Returns:
          Dictionary of float horizontal positions, keyed by node id.
        """
        # Statistics needed for horizontal positions.
        number_descendants = dict.fromkeys(self.get_node_list(), 1)
        # number_descendants includes the key node itself
        horizontal_lower_bound = dict.fromkeys(self.get_node_list(), 0.0)
        horizontal_upper_bound = dict.fromkeys(self.get_node_list(), 1.0)
        horizontal_positions = dict.fromkeys(self.get_node_list())
        visited = dict.fromkeys(self.get_node_list(), False)
        # Count the number of descendants for each node.
        # Do a post-order traversal of the tree.
        node_stack = []
        node_stack.append(self.root.name)
        while node_stack:
            current_node = node_stack[len(node_stack) - 1]
            lchild = self.get_left_child(current_node)
            rchild = self.get_right_child(current_node)
            is_node_added = False
            # Add the next unvisited child to the stack
            if lchild is not None and not visited[quote(lchild)]:
                node_stack.append(lchild)
                is_node_added = True
            if (rchild is not None and not visited[quote(rchild)] and
                is_node_added==False):
                node_stack.append(rchild)
                is_node_added = True
            # If all childs visited, then update number_descendants
            if not is_node_added:
                if lchild is not None:
                    number_descendants[quote(current_node)] += (
                                number_descendants[quote(lchild)])
                if rchild is not None:
                    number_descendants[quote(current_node)] += (
                                number_descendants[quote(rchild)])
                visited[quote(current_node)] = True
                del node_stack[len(node_stack) - 1]
        # Traverse the tree and set horizontal positions.
        # Do a pre-order traversal of the tree.
        node_stack = []
        node_stack.append(self.root.name)
        horizontal_lower_bound[self.root.name] = 0.0
        horizontal_upper_bound[self.root.name] = 1.0
        while node_stack:
            node = node_stack.pop()
            lchild = self.get_left_child(node)
            rchild = self.get_right_child(node)
            direction = None
            number_of_children = 0
            children_list = []
            # Place all children on the stack
            if lchild is not None:
                node_stack.append(lchild)
                number_of_children += 1
                direction = 'L'
                children_list.append(lchild)
            if rchild is not None:
                node_stack.append(rchild)
                number_of_children += 1
                direction = 'R'
                children_list.append(rchild)
            # Convenience variables
            current_lower_bound = horizontal_lower_bound[quote(node)]
            current_range = (horizontal_upper_bound[quote(node)] -
                             horizontal_lower_bound[quote(node)])
            total_descendants = number_descendants[quote(node)]
            sorted_child_labels = sorted(children_list)
            # Determine where to place this node with respect to its children.
            # Put the node in the center, or have more children on the left.
            before_index = int(math.ceil(old_div(number_of_children,2.0)))
            # Exception with a single node that is 'L'
            if number_of_children == 1:
                if direction != 'L':
                    before_index = 0
            cumulative_descendants = 0
            for i, label in enumerate(sorted_child_labels):
                if before_index == i:
                    # Determine the relative position for the current node
                    relative_position = old_div((cumulative_descendants + 0.5), (
                        total_descendants))
                    cumulative_descendants += 1
                # Set bounds for this child
                horizontal_lower_bound[quote(label)] = (
                    current_lower_bound + float(cumulative_descendants) /
                    total_descendants * current_range)
                # Increment cumulative_descendants, which also lets us compute
                # the upper bound.
                cumulative_descendants += number_descendants[quote(label)]
                horizontal_upper_bound[quote(label)] = (
                    current_lower_bound + float(cumulative_descendants) /
                    total_descendants * current_range)
            # Catch the case that the node comes after all its children.
            # This must also work for the case that this is the only node.
            if before_index == len(sorted_child_labels):
                relative_position = old_div((cumulative_descendants + 0.5), (
                    total_descendants))
            # Finally set the position for the current node
            horizontal_positions[quote(node)] = (
                horizontal_lower_bound[quote(node)] + relative_position * (
                    horizontal_upper_bound[quote(node)] -
                    horizontal_lower_bound[quote(node)]))
        return horizontal_positions

    def WriteDataFileFromList(self, filename, data_list):
        """
        Write a list of string data to a file with one entry per line.
        Args:
          filename: String filename to open.
          data_list: List of string values to write.
        """
        outfile = open(filename, 'w')
        for line in data_list:
            outfile.write(line)
        outfile.close()

    def GenerateTreeImage(self, fixed_horizontal_positions = False):
        """
        Generate files necessary for a tree image.
        Two files are necessary: a data file and a Gnuplot script file (which
        references the data file).
        """
        index_string = self.GetImageCounterString()
        if fixed_horizontal_positions:
            name_prefix = 'tree_alt'
            horizontal_positions = self.GetTreeFixedHorizontalPositions()
        else:
            name_prefix = 'tree'
            horizontal_positions = self.GetTreeHorizontalPositions()
        candidate_lines = []
        pregnant_lines = []
        branched_lines = []
        infeasible_lines = []
        fathomed_lines = []
        integer_lines = []
        additional_script_lines = []
        node_list = self.get_node_list()
        print_edges = (len(node_list) <= self._edge_limit)
        for node in node_list:
            node_lp_bound = self.get_node_attr(node, 'lp_bound')
            if self.get_node_attr(node, 'status') == 'candidate':
                # TODO(bhunsaker): add optional fathoming check
                candidate_lines.append('%0.6f %0.6f\n' % (
                        horizontal_positions[node], node_lp_bound))
            elif self.get_node_attr(node, 'status') == 'pregnant':
                # TODO(bhunsaker): add optional fathoming check
                pregnant_lines.append('%0.6f %0.6f\n' % (
                        horizontal_positions[node], node_lp_bound))
            elif self.get_node_attr(node, 'status') == 'branched':
                branched_lines.append('%0.6f %0.6f\n' % (
                        horizontal_positions[node], node_lp_bound))
            elif self.get_node_attr(node, 'status') == 'infeasible':
                infeasible_lines.append('%0.6f %0.6f\n' % (
                        horizontal_positions[node], node_lp_bound))
            elif self.get_node_attr(node, 'status') == 'fathomed':
                fathomed_lines.append('%0.6f %0.6f\n' % (
                        horizontal_positions[node], node_lp_bound))
            elif self.get_node_attr(node, 'status') == 'integer':
                integer_lines.append('%0.6f %0.6f\n' % (
                        horizontal_positions[node], node_lp_bound))
            if print_edges and node != self.root.name:
                if True:
                    _parent_id = self.get_node_attr(node, 'parent')
                    additional_script_lines.append(
                     'set arrow from %0.6f, %0.6f to %0.6f, %0.6f nohead lt -1 '
                     'lw 0.2\n' % (horizontal_positions[_parent_id],
                     self.get_node_attr(_parent_id, 'lp_bound'),
                     horizontal_positions[node],
                     self.get_node_attr(node, 'lp_bound')))
        plot_parts = []
        # Plot root node.
        plot_parts.append('"< echo %0.6f %0.6f" w p lt 2 pt 7' %
                          (horizontal_positions[self.root.name],
                           self.root.get_attr('lp_bound')))
        # If desired, sample from the set of nodes rather than plotting all.
        if self._sample_tree:
            sample_size = self._sample_tree
            if len(branched_lines) > sample_size:
                branched_lines = random.sample(branched_lines, sample_size)
            if len(fathomed_lines) > sample_size:
                fathomed_lines = random.sample(fathomed_lines, sample_size)
            if len(infeasible_lines) > sample_size:
                infeasible_lines = random.sample(infeasible_lines, sample_size)
            if len(pregnant_lines) > sample_size:
                pregnant_lines = random.sample(pregnant_lines, sample_size)
            if len(candidate_lines) > sample_size:
                candidate_lines = random.sample(candidate_lines, sample_size)
            if len(integer_lines) > sample_size:
                integer_lines = random.sample(integer_lines, sample_size)
        # Output all data files.  Note that the order below matters.
        if len(branched_lines):
            self.WriteDataFileFromList('%s_branched%s.dat' % (name_prefix,
                                                              index_string),
                                       branched_lines)
            plot_parts.append('\'%s_branched%s.dat\' w p lt rgb "yellow" pt 7' %
                              (name_prefix, index_string))
        if len(fathomed_lines):
            self.WriteDataFileFromList('%s_fathomed%s.dat' % (name_prefix,
                                                              index_string),
                                       fathomed_lines)
            plot_parts.append('\'%s_fathomed%s.dat\' w p lt rgb "light-red" pt 7' %
                              (name_prefix, index_string))
        if len(infeasible_lines):
            self.WriteDataFileFromList('%s_infeasible%s.dat' % (name_prefix,
                                                                index_string),
                                       infeasible_lines)
            plot_parts.append('\'%s_infeasible%s.dat\' w p lt rgb "dark-red" pt 7' %
                              (name_prefix, index_string))
        if len(pregnant_lines):
            self.WriteDataFileFromList('%s_pregnant%s.dat' % (name_prefix,
                                                              index_string),
                                       pregnant_lines)
            plot_parts.append('\'%s_pregnant%s.dat\' w p lt rgb "green" pt 7' %
                              (name_prefix, index_string))
        if len(candidate_lines):
            for line in candidate_lines:
                plot_parts.append('"< echo %s" w p lt rgb "green" pt 7'
                                  %line.rstrip('\r\n'))
        if len(integer_lines):
            self.WriteDataFileFromList('%s_integer%s.dat' % (name_prefix,
                                                             index_string),
                                       integer_lines)
            plot_parts.append('\'%s_integer%s.dat\' w p lt rgb "cyan" pt 7' %
                              (name_prefix, index_string))
        if self._incumbent_value is not None:
            plot_parts.append('%0.6f lt 1 lw 0.5' % self._incumbent_value)
        additional_script_lines.append('plot %s\n' % ', '.join(plot_parts))
        additional_script_lines.append('unset arrow\n')
        image_min_obj, image_max_obj = self.GetImageObjectiveBounds(
            self._min_objective_value, self._max_objective_value)
        data = ''
        data += 'set terminal png notransparent size 480,360\n'
        data += 'set nokey\n'
        data += 'set autoscale\n'
        data += 'set tics scale 0.001\n'
        data += 'set pointsize 0.5\n'
        data += 'set xrange [-0.1:1.1]\n'
        data += 'set yrange [%0.6f:%0.6f]\n' % (image_max_obj,
                                                          image_min_obj)
        data += 'set format x ""\n'
        data += 'set ylabel "obj. value"\n'
        if self._filename is None:
            data += 'set title "B&B tree"\n\n'
        else:
            data += 'set title "B&B tree (%s %.2fs %s)"\n\n' % (
                self._filename, self._time, self._label)
        for line in additional_script_lines:
            data += line

        return data

    def ProcessLine(self, line):
        """
        Process a line of the input file, generating images if appropriate.
        Parses the line, updates internal data structures, and creates images
        if appropriate.
        Args:
          line: String input line to process.
        """
        line = line.strip()
        # Comments start with a '#'
        if line[0] == '#':
            return
        tokens = line.split()
        if len(tokens) < 3:
            print('Incomplete or invalid line: %s' %' '.join(tokens))
            sys.exit(1)
        # Tokens shared by all line types
        self._time = float(tokens[0])
        line_type = tokens[1]
        remaining_tokens = tokens[2:]
        # Process the line based on the type
        if line_type == 'heuristic':
            self._optimal_soln_time = self._time
            self.ProcessHeuristicLine(remaining_tokens)
        else:
            # Other node types share common tokens
            node_id = int(tokens[2])
            parent_id = int(tokens[3])
            branch_direction = tokens[4]
            remaining_tokens = tokens[5:]
            # TODO(aykut):parent id of root node is 0 when we read from file.
            if id==self.root:
                parent_id = None
            # Check that the parent node id is valid
            # elif parent_id not in self.get_node_list() and self.root is not None:
            #     print 'Parent id does not exist: %s' % line
            #     sys.exit(1)
            if line_type == 'integer':
                self._optimal_soln_time = self._time
                self.ProcessIntegerLine(node_id, parent_id,
                                        branch_direction, remaining_tokens)
            elif line_type == 'fathomed':
                self.ProcessFathomedLine(node_id, parent_id,
                                         branch_direction, remaining_tokens)
            elif line_type == 'candidate':
                self.ProcessCandidateLine(node_id, parent_id,
                                          branch_direction, remaining_tokens)
            elif line_type == 'pregnant':
                self.ProcessPregnantLine(node_id, parent_id,
                                         branch_direction, remaining_tokens)
            elif line_type == 'branched':
                self.ProcessBranchedLine(node_id, parent_id,
                                         branch_direction, remaining_tokens)
            elif line_type == 'infeasible':
                self.ProcessInfeasibleLine(node_id, parent_id,
                                           branch_direction, remaining_tokens)
            else:
                print('Unexpected line type "%s": %s' % (line_type,
                                                         ' '.join(tokens)))
                sys.exit(1)

    def ProcessHeuristicLine(self, remaining_tokens):
        """
        Core processing for a line of type 'heuristic'.
        Args:
          remaining_tokens: List of string tokens. These are those that remain
            after any common tokens are processed.
        """
        # Parse remaining tokens
        if len(remaining_tokens) < 1 or len(remaining_tokens) > 2:
            print('Invalid line: %s heuristic %s' % (
                    self._time, ' '.join(remaining_tokens)))
            print('Should match: <time> heuristic <obj value>'+\
                ' [<associated node id>]')
            sys.exit(1)
        objective_value = float(remaining_tokens[0])
        if len(remaining_tokens) == 2:
            associated_node = remaining_tokens[1]
        else:
            associated_node = None
        # Check that this is actually an improvement
        if self._incumbent_value is not None and self._optimization_sense is None:
            if objective_value > self._incumbent_value:
                print("Objective sense unset, guessing maximization")
                self._optimization_sense = 'max'
            else:
                print("Objective sense unset, guessing minimization")
                self._optimization_sense = 'min'
        if not self.IsBetterThanIncumbent(objective_value):
            return
        self._previous_incumbent_value = self._incumbent_value
        self._incumbent_value = objective_value
        self.UpdateObjectiveValueLimits(objective_value)
        self._incumbent_parent = associated_node
        # Set variable to generate images
        self._new_integer_solution = True

    def ProcessIntegerLine(self, node_id, parent_id, branch_direction,
                           remaining_tokens):
        """
        Core processing for a line of type 'integer'.
        Args:
          node_id: String node id.
          parent_id: String node id of parent.
          branch_direction: String of 'L' or 'R' indicating whether this node
          is the left or right child of its parent.
          remaining_tokens: List of string tokens. These are those that remain
            after any common tokens are processed.
        """
        # Parse remaining tokens
        if len(remaining_tokens) != 1:
            print('Invalid line: %s integer %s %s %s %s' % (
                    self._time, node_id, parent_id, branch_direction,
                    ' '.join(remaining_tokens)))
            print('Should match: <time> integer <node id> <parent id>'+\
                '<branch direction> <obj value>')
            sys.exit(1)
        objective_value = float(remaining_tokens[0])
        self.AddOrUpdateNode(node_id, parent_id, branch_direction, 'integer',
                             objective_value, None, None)
        self._previous_incumbent_value = self._incumbent_value
        self._incumbent_value = objective_value
        self._incumbent_parent = parent_id
        self._new_integer_solution = True

    def ProcessFathomedLine(self, node_id, parent_id, branch_direction,
                            remaining_tokens):
        """
        Core processing for a line of type 'fathomed'.
        Args:
          node_id: String node id.
          parent_id: String node id of parent.
          branch_direction: String of 'L' or 'R' indicating whether this node is
            the left or right child of its parent.
          remaining_tokens: List of string tokens. These are those that remain
            after any common tokens are processed.
        """
        # Print a warning if there is no current incumbent.
        if self._incumbent_value is None:
            print('WARNING: Encountered "fathom" line before first incumbent.')
            print('  This may indicate an error in the input file.')
        # Parse remaining tokens
        if len(remaining_tokens) > 1:
            print('Invalid line: %s fathomed %s %s %s %s' % (
                    self._time, node_id, parent_id, branch_direction,
                    ' '.join(remaining_tokens)))
            print('Should match: <time> fathomed <node id> <parent id>'+\
                '<branch direction> [<lp bound>]')
            sys.exit(1)
        if len(remaining_tokens) == 1:
            lp_bound = float(remaining_tokens[0])
        else:
            if (node_id in self.get_node_list() and
                self.get_node_attr(node_id, 'lp_bound') is not None):
                lp_bound = self.get_node_attr(node_id, 'lp_bound')
            else:
                lp_bound = self.get_node_attr(parent_id, 'lp_bound')
            if self._optimization_sense == 'min':
                if (self._incumbent_value is not None and
                    lp_bound < self._incumbent_value):
                    lp_bound = self._incumbent_value
            elif self._optimization_sense == 'max':
                if (self._incumbent_value is not None and
                    lp_bound > self._incumbent_value):
                    lp_bound = self._incumbent_value
        parent_node = self.get_node(parent_id)
        self.AddOrUpdateNode(node_id, parent_id, branch_direction, 'fathomed',
                             lp_bound,
                             self.get_node_attr(parent_id,
                                                'integer_infeasibility_count'),
                             self.get_node_attr(parent_id,
                                                'integer_infeasibility_sum'))

    def ProcessPregnantLine(self, node_id, parent_id, branch_direction,
                            remaining_tokens):
        """
        Core processing for a line of type 'pregnant'.
        Args:
          node_id: String node id.
          parent_id: String node id of parent.
          branch_direction: String of 'L' or 'R' indicating whether this node is
            the left or right child of its parent.
          remaining_tokens: List of string tokens. These are those that remain
            after any common tokens are processed.
        """
        # Parse remaining tokens
        if len(remaining_tokens) != 3:
            print('Invalid line: %s pregnant %s %s %s %s' % (
                    self._time, node_id, parent_id, branch_direction,
                    ' '.join(remaining_tokens)))
            print('Should match: <time> pregnant <node id> <parent id> ')
            print('<branch direction> <lp bound> ')
            print('<sum of integer infeasibilities> <number of integer ')
            print('infeasibilities>')
            sys.exit(1)
        lp_bound = float(remaining_tokens[0])
        integer_infeasibility_sum = float(remaining_tokens[1])
        integer_infeasibility_count = int(remaining_tokens[2])

        self.AddOrUpdateNode(node_id, parent_id, branch_direction, 'pregnant',
                             lp_bound, integer_infeasibility_count,
                             integer_infeasibility_sum)

    def ProcessBranchedLine(self, node_id, parent_id, branch_direction,
                            remaining_tokens):
        """
        Core processing for a line of type 'branched'.
        Args:
          node_id: String node id.
          parent_id: String node id of parent.
          branch_direction: String of 'L' or 'R' indicating whether this node
          is the left or right child of its parent.
          remaining_tokens: List of string tokens. These are those that remain
            after any common tokens are processed.
        """
        # Parse remaining tokens
        if len(remaining_tokens) not in [3, 5]:
            print('Invalid line: %s branched %s %s %s %s' % (
                    self._time, node_id, parent_id, branch_direction,
                    ' '.join(remaining_tokens)))
            print('Should match: <time> branched <node id> <parent id> ')
            print('<branch direction> <lp bound> ')
            print('<sum of integer infeasibilities> <number of integer ')
            print('infeasibilities>')
            sys.exit(1)
        lp_bound = float(remaining_tokens[0])
        integer_infeasibility_sum = float(remaining_tokens[1])
        integer_infeasibility_count = int(remaining_tokens[2])
        condition_begin = None
        condition_end = None
        if len(remaining_tokens) == 5:
            # In this case, we must also be printing conditions numbers
            condition_begin = int(remaining_tokens[3])
            condition_end = int(remaining_tokens[4])
        self.AddOrUpdateNode(node_id, parent_id, branch_direction, 'branched',
                             lp_bound, integer_infeasibility_count,
                             integer_infeasibility_sum, condition_begin,
                             condition_end)

    def ProcessInfeasibleLine(self, node_id, parent_id, branch_direction,
                              remaining_tokens):
        """
        Core processing for a line of type 'infeasible'.
        Args:
          node_id: String node id.
          parent_id: String node id of parent.
          branch_direction: String of 'L' or 'R' indicating whether this node is
            the left or right child of its parent.
          remaining_tokens: List of string tokens. These are those that remain
            after any common tokens are processed.
        """
        # Parse remaining tokens
        if len(remaining_tokens) not in [0, 2]:
            print('Invalid line: %s infeasible %s %s %s %s' % (
                    self._time, node_id, parent_id, branch_direction,
                    ' '.join(remaining_tokens)))
            print('Should match: <time> infeasible <node id> <parent id> ')
            print('<branch direction>')
            sys.exit(1)
        # Use parent values if the node does not have its own
        lp_bound = self.get_node_attr(parent_id, 'lp_bound')
        ii_count = self.get_node_attr(parent_id, 'integer_infeasibility_count')
        ii_sum = self.get_node_attr(parent_id, 'integer_infeasibility_sum')
        if node_id in self.get_node_list():
            if self.get_node_attr(node_id, 'lp_bound') is not None:
                lp_bound = self.get_node_attr(node_id, 'lp_bound')
            if (self.get_node_attr(node_id, 'integer_infeasibility_count')
                is not None):
                ii_count = self.get_node_attr(node_id,
                                              'integer_infeasibility_count')
            if (self.get_node_attr(node_id, 'integer_infeasibility_sum')
                is not None):
                ii_sum = self.get_node_attr(node_id,'integer_infeasibility_sum')
        if len(remaining_tokens) == 2:
            # In this case, we must also be printing conditions numbers
            condition_begin = int(remaining_tokens[0])
            condition_end = int(remaining_tokens[1])
        self.AddOrUpdateNode(node_id, parent_id, branch_direction, 'infeasible',
                             lp_bound, ii_count, ii_sum)

    def ProcessCandidateLine(self, node_id, parent_id, branch_direction,
                             remaining_tokens):
        """
        Core processing for a line of type 'candidate'.
        Args:
          node_id: String node id.
          parent_id: String node id of parent.
          branch_direction: String of 'L' or 'R' indicating whether this node
          is the left or right child of its parent.
          remaining_tokens: List of string tokens. These are those that remain
            after any common tokens are processed.
        """
        # Parse remaining tokens
        if len(remaining_tokens) == 2 or len(remaining_tokens) > 3:
            print('Invalid line: %s branched %s %s %s %s' % (
                    self._time, node_id, parent_id, branch_direction,
                    ' '.join(remaining_tokens)))
            print('Should match: <time> candidate <node id> <parent id> ')
            print('<branch direction> [<lp bound>] ')
            print('[<sum of integer infeasibilities> <number of integer ')
            print('infeasibilities>]')
            sys.exit(1)
        # if parent_id not in self.get_node_list():
        #     print 'Error: node %s not in set' % parent_id
        #     sys.exit(1)
        # TODO(bhunsaker): Check that we handle the cases of updating a
        #candidate.
        if len(remaining_tokens) > 0:
            lp_bound = float(remaining_tokens[0])
        else:
            lp_bound = self.get_node_attr(parent_id, 'lp_bound')
        if len(remaining_tokens) == 3:
            integer_infeasibility_sum = float(remaining_tokens[1])
            integer_infeasibility_count = int(remaining_tokens[2])
        else:
            integer_infeasibility_sum = self.get_node_attr(parent_id,
                                                  'integer_infeasibility_sum')
            integer_infeasibility_count = self.get_node_attr(parent_id,
                                                'integer_infeasibility_count')
        self.AddOrUpdateNode(node_id, parent_id, branch_direction, 'candidate',
                             lp_bound, integer_infeasibility_count,
                             integer_infeasibility_sum)

    def RunGnuplotOnAllFiles(self):
        """Runs Gnuplot on all files in self._gnuplot_files."""
        for file in self._gnuplot_files:
            subprocess.call(['gnuplot', file])

    def CreateAnimatedImages(self):
        """Create animated images based on the static images."""
        histogram_re = re.compile('histogram')
        histogram_images = [re.sub('gnuplot', 'png', file)
                            for file in self._gnuplot_files
                            if histogram_re.match(file)]
        if len(histogram_images):
            args = ['convert', '-delay', '15', '-loop', '1']
            args.extend(histogram_images)
            args.append('animated_histogram.gif')
            subprocess.call(args)
        scatterplot_re = re.compile('scatterplot')
        scatterplot_images = [re.sub('gnuplot', 'png', file)
                              for file in self._gnuplot_files
                              if scatterplot_re.match(file)]
        if len(scatterplot_images):
            args = ['convert', '-delay', '15', '-loop', '1']
            args.extend(scatterplot_images)
            args.append('animated_scatterplot.gif')
            subprocess.call(args)
        tree_re = re.compile('tree\.')
        tree_images = [re.sub('gnuplot', 'png', file)
                       for file in self._gnuplot_files
                       if tree_re.match(file)]
        if len(tree_images):
            args = ['convert', '-delay', '15', '-loop', '1']
            args.extend(tree_images)
            args.append('animated_tree.gif')
            subprocess.call(args)
        tree_alt_re = re.compile('tree_alt')
        tree_alt_images = [re.sub('gnuplot', 'png', file)
                           for file in self._gnuplot_files
                           if tree_alt_re.match(file)]
        if len(tree_alt_images):
            args = ['convert', '-delay', '15', '-loop', '1']
            args.extend(tree_alt_images)
            args.append('animated_tree_alt.gif')
            subprocess.call(args)

    def GeneratePredictionImages(self):
        gap_measures = self._objective_gap_forecaster.GetAllMeasures()
        ssg_measures = self._sum_subtree_gaps_forecaster.GetAllMeasures()
        # Check that there are values to process.
        if len(gap_measures) == 0 or len(ssg_measures) == 0:
            print('WARNING: Not printing prediction images because at least'+\
                ' one measure set is empty.')
            print('  Gap measures: %d' % len(gap_measures))
            print('  SSG measures: %d' % len(ssg_measures))
            return
        # Gap measures
        gap_data_filename = 'gap_measures.dat'
        data_file = open(gap_data_filename, 'w')
        for measure in gap_measures:
            data_file.write('%0.6f %0.6f\n' % (measure.time, measure.value))
        data_file.close()
        # SSG measures
        ssg_data_filename = 'ssg_measures.dat'
        data_file = open(ssg_data_filename, 'w')
        # We need to scale the SSG measures so that it will make sense to
        # look at them on the same plot with gap measures.
        scale_factor=old_div(float(gap_measures[0].value),float(ssg_measures[0].value))
        for measure in ssg_measures:
            data_file.write('%0.6f %0.6f\n' % (measure.time,
                                               measure.value * scale_factor))
        data_file.close()
        # Set terminal for the output files.
        measures_script = 'set terminal png notransparent size 480,360\n\n'
        # Make settings for the plot.
        if self.filename is None:
            measures_script += 'set title "Progress Measures"\n'
        else:
            measures_script += ('set title "Progress Measures: %s, %s"\n' % (
                    self._filename, self._label))
        measures_script += 'set xlabel \"time (s)\"\n'
        measures_script += 'set ylabel \"measure\"\n'
        measures_script += 'set autoscale\n'
        # Plot the data points.
        measures_script += (
            'plot \'%s\' with linespoints linetype 3 title \"(SSG)\", '
            '\'%s\' with linespoints linetype 4 pointtype 19 '
            'title \"(MIP gap)\"\n' %
            (ssg_data_filename, gap_data_filename))
        measures_script += 'show output\n'

        return measures_script

    def GenerateForecastImages(self):
        # Forecasts
        # Gap forecasts
        gap_forecasts = self._objective_gap_forecaster.GetAllForecasts()
        gap_data_filename = 'gap_forecasts.dat'
        if gap_forecasts:
            data_file = open(gap_data_filename, 'w')
            for forecast in gap_forecasts:
                data_file.write('%0.6f %0.6f\n' % (forecast.time,
                                                   forecast.forecast))
            data_file.close()
        # SSG forecasts
        ssg_forecasts = self._sum_subtree_gaps_forecaster.GetAllForecasts()
        ssg_data_filename = 'ssg_forecasts.dat'
        if ssg_forecasts:
            data_file = open(ssg_data_filename, 'w')
            for forecast in ssg_forecasts:
                data_file.write('%0.6f %0.6f\n' % (forecast.time,
                                                   forecast.forecast))
            data_file.close()
        if not gap_forecasts and not ssg_forecasts:
            print('No forecasts made, so not creating forecast image.')
            return
        # Set terminal for the output files.
        forecast_script = 'set terminal png notransparent size 480,360\n\n'
        # Make settings for the plot.
        if self._filename is None:
            forecast_script += 'set title "Forecasts"\n'
        else:
            forecast_script += ('set title "Forecasts: %s, %s"\n' % (
                    self._filename, self._label))
        forecast_script += 'set xlabel \"time (s)\"\n'
        forecast_script += 'set ylabel \"prediction of total time\"\n'
        forecast_script += 'set autoscale\n'
        # Plot the data points and the unit-slope line (to show elapsed time).
        forecast_script += 'plot '
        if forecast_forecasts:
            forecast_script += ('\'%s\' with linespoints linetype 3 '
                              'title \"(SSG)\", ' % ssg_data_filename)
        if gap_forecasts:
            forecast_script += ('\'%s\' with linespoints linetype 4 pointtype 19 '
            'title \"(MIP gap)\", ' % gap_data_filename)
        forecast_script += 'x linetype 0 title \"elapsed time\"\n'
        forecast_script += 'show output\n'

        return forecast_script

    def _get_fh(self, path, mode='r'):
        '''
        Return a file handle for given path.
        Path can be a string or a file handle.
        Attempt to uncompress/compress files ending in '.gz' and '.bz2'.
        '''
        if self._is_string_like(path):
            if path.endswith('.gz'):
#               import gzip
#               fh = gzip.open(path,mode=mode)  # doesn't return real fh
                fh=os.popen("gzcat "+path) # probably not portable
            elif path.endswith('.bz2'):
#               import bz2
#               fh = bz2.BZ2File(path,mode=mode) # doesn't return real fh
                fh=os.popen("bzcat "+path) # probably not portable
            else:
                fh = file(path,mode=mode)
        elif hasattr(path, 'write'):
            # Note, mode of file handle is unchanged.
            fh = path
        else:
            raise TypeError('path must be a string or file handle.')
        return fh

    def _is_string_like(self, obj): # from John Hunter, types-free version
        try:
            obj + ''
        except (TypeError, ValueError):
            return False
        return True

def CreatePerlStyleBooleanFlag(parser, flag_text, variable_name, help_text):
    """
    Add two options to an optparse.OptionParser, one with a 'no' prefix.
    Two options are created.  One has the flag_text and one has 'no' prepended
    to the flag_text.  For example, --foo and --nofoo.  This is similar to a
    common style in Perl.
    Args:
      parser: optparse.OptionParser object.
      flag_text: String text for the flag.
      variable_name: String name of the variable to store the flag results.
      help_text: String that describes the flag.
    """
    parser.add_option('--' + flag_text,
                      action='store_true', dest=variable_name, default=False,
                      help=help_text)
    parser.add_option('--no' + flag_text,
                      action='store_false', dest=variable_name,
                      help='do not ' + flag_text)

def parse_options():
    ''' Parse arguments and flags'''
    usage_text = 'usage: %prog [options] <input file>'
    parser = optparse.OptionParser(usage=usage_text)

    parser.add_option('--interval', dest='interval', type='float',
                      default=10.0,
                      help='generate images every TIME seconds',
                      metavar='TIME')
    parser.add_option('--time_limit', dest='time_limit', type='float',
                      help='process at most TIME seconds of solver time',
                      metavar='TIME')
    parser.add_option('--label', dest='label', default='',
                      help='add LABEL to all images',
                      metavar='LABEL')
    parser.add_option('--use_common_bounds', dest='use_common_bounds',
                      action='store_true', default=False,
                      help='use the same bounds on all images')
    CreatePerlStyleBooleanFlag(parser, 'fathom', 'fathom',
                               'do an interval check for fathomed nodes')
    CreatePerlStyleBooleanFlag(parser, 'histogram', 'histogram',
                               'print histograms')
    CreatePerlStyleBooleanFlag(parser, 'scatterplot', 'scatterplot',
                               'print scatterplots')
    CreatePerlStyleBooleanFlag(parser, 'path', 'path',
                               'print scatterplot incumbent paths')
    CreatePerlStyleBooleanFlag(parser, 'tree', 'tree',
                               'print tree images')
    CreatePerlStyleBooleanFlag(parser, 'fixed_tree', 'fixed_tree',
                               'print tree images with fixed horizontal '
                               'positions')
    CreatePerlStyleBooleanFlag(parser, 'predictions', 'predictions',
                               'print time predictions')
    parser.add_option('--all', dest='all_images', action='store_true',
                      default=False, help='print all images')
    parser.add_option('--logscaley', dest='logscaley',
                      action='store_true', default=False,
                      help='use log scale for histogram sizes')
    parser.add_option('--animate', dest='animate',
                      action='store_true', default=False,
                      help='create animated GIF of each image set')
    parser.add_option('--no_run_gnuplot', dest='run_gnuplot',
                      action='store_false', default=True,
                      help='do not run Gnuplot on the generated files')
    parser.add_option('--edge_limit', dest='edge_limit', type='int',
                      default=1000000,
                      help='do not print edges in tree plots if more than '
                      'NUM nodes',
                      metavar='NUM')
    parser.add_option('--sample_tree', dest='sample_tree', type='int',
                      default=0,
                      help='use at most NUM nodes of each type in tree images;'
                      ' zero means no limit',
                      metavar='NUM')
    (options, args) = parser.parse_args()
    if (len(args) != 1):
        parser.print_usage()
        sys.exit(1)
    input_filename = args[0]
    if options.all_images:
        options.histogram = True
        options.scatterplot = True
        options.path = True
        options.tree = True
        options.fixed_tree = True
        options.predictions = True
    # Abort if no images chosen
    if (not options.histogram and not options.scatterplot and
        not options.path and not options.tree and not options.fixed_tree and
        not options.predictions):
        print('No image types specified so not processing.')
        parser.print_usage()
        sys.exit(1)
    # Bounds for incumbent paths will be undefined without scatterplots.
    if options.path:
        options.scatterplot = True
    # TODO(bhunsaker): Check whether Gnuplot can be accessed early.
    # TODO(bhunsaker): Check whether animation program can be accessed.
    if (len(args) != 1):
        parser.print_usage()
        sys.exit(1)
    if options.all_images:
        options.histogram = True
        options.scatterplot = True
        options.path = True
        options.tree = True
        options.fixed_tree = True
        options.predictions = True
    # Abort if no images chosen
    if (not options.histogram and not options.scatterplot and
        not options.path and not options.tree and not options.fixed_tree and
        not options.predictions):
        print('No image types specified so not processing.')
        sys.exit(1)
    # Bounds for incumbent paths will be undefined without scatterplots.
    if options.path:
        options.scatterplot = True
    return (input_filename, options)


if __name__ == '__main__':
    from .BranchAndBound import GenerateRandomMIP, BranchAndBound

    T = BBTree()
    #T.set_layout('dot2tex')
    #T.set_display_mode('file')
    T.set_display_mode('matplotlib')
    CONSTRAINTS, VARIABLES, OBJ, MAT, RHS = GenerateRandomMIP(rand_seed = 120)
    BranchAndBound(T, CONSTRAINTS, VARIABLES, OBJ, MAT, RHS,
                   branch_strategy = MOST_FRACTIONAL,
                   search_strategy = BEST_FIRST,
                   display_interval = 10000)
