from future import standard_library
standard_library.install_aliases()
from builtins import str
try:
    from coinor.grumpy import BBTree
except ImportError:
    from src.grumpy import BBTree
import sys
import io
from PIL import Image as PIL_Image

bt = BBTree()
bt.set_display_mode('pygame')
line_number = 1
file_ = open('p0201_GLPK.vbc', 'r')
#for line in sys.stdin:
for line in file_:
    bt.ProcessLine(line)
    #To print out snapshots of the tree
    if line_number%100000 != 0:
        continue
    imagefile = open('tree-'+str(line_number)+'.png','w')
    imagefile.write(bt.GenerateTreeImage())
    imagefile.close()

#gnuplot_image = bt.GenerateHistogram()
gnuplot_image = io.StringIO(bt.GenerateTreeImage())
#gnuplot_image = bt.GenerateScatterplot()
#gnuplot_image = bt.GenerateIncumbentPath()
#gnuplot_image = bt.GenerateForecastImages()
im = PIL_Image.open(gnuplot_image)
im.show()
#bt.display_all()

bt.set_display_mode('file')
bt.display()
