#GrUMPy 0.8

Graphics for Understanding Mathematical Programming in Python (GrUMPy) is a
Python library for visualizing various aspects of mathematical programming,
including visualizations of the branch-and process, branch-and-bound trees,
polyhedra, cutting plane methods, etc. The goal is clarity in implementation
rather than eficiency. Most methods have an accompanying visualization and are
thus appropriate for use in the classroom.

Documentation for the API is here:

http://pythonhosted.org/coinor.grumpy

Installation:

easy_install coinor.grumpy

##Installation Notes

1. GrUMPy depends on [GiMPy](https://github.com/coin-or/GiMPy), which will be 
  automatically installed as part of the setup. However, in order for GiMPy to
  visualize the branch-and-bound tree, it's necessary to install 
  [GraphViz](http://www.graphviz.org/Download.php) and choose one of these 
  additional methods for display:
  * Recommanded: [xdot](https://pypi.python.org/pypi/xdot) along with 
    [PyGtk](http://www.pygtk.org/) and call `set_display_mode('xdot')`
  * [Python Imaging Library](http://www.pythonware.com/products/pil/) and 
    call `set_display_mode('PIL')`
  * [Pygame](pygame.org) and call `set_display_mode('pygame')`
  * Call `set_display_mode('file')` to just write files to disk that have to
    then be opened manually. 
  
  It is also possible to typeset labels in LaTex and to output the graph in 
  LaTex format using `dot2tex`. After installing `dot2tex`, this can be done 
  by simply calling the method `write(basename='fileName', format='dot')`, and 
  then doing `dot2tex --tmath fileName.dot` or by calling 
  `set_display_mode('dot2tex')` and then `display()` as usual. At the moment,
  the latter only seems to work with version `2.9.0dev` available 
  [here](https://github.com/Alwnikrotikz/dot2tex). For the former method, just 
  using `easy_install dot2tex` should work fine.
2. GrUMPy can also visualize 2D polyhedra with the installation of 
  [pypolyhedron](https://github.com/rdeits/pypolyhedron), which must be
  installed from source.

##Examples

![Branch and bound tree](https://github.com/coin-or/GrUMPy/raw/master/images/BranchAndBound.png)

![Figure 1 made with GrUMPy](https://raw.githubusercontent.com/tkralphs/GrUMPy/master/images/polyhedron.png)

![Figure 2 made with GrUMPy](https://raw.githubusercontent.com/tkralphs/GrUMPy/master/images/GMI-Row2-Disjunction.png)

![Figure 3 made with GrUMPy](https://raw.githubusercontent.com/tkralphs/GrUMPy/master/images/GMI-Row3.png)
