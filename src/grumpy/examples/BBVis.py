from future import standard_library
standard_library.install_aliases()
from builtins import str
import coinor.grumpy.examples
try:
    from coinor.grumpy import BBTree
except ImportError:
    from src.grumpy import BBTree
from inspect import getfile
from os.path import join, dirname
import sys
import io

bt = BBTree()
line_number = 0
#for line in sys.stdin:
instance = join(dirname(getfile(coinor.grumpy.examples)),
                'p0201_GLPK.vbc')
with open(instance, 'r') as file_:
    bt.set_display_mode('file')
    for line in file_:
        line_number += 1
        bt.ProcessLine(line)
        #To print out snapshots of the tree
        if line_number%100 != 0:
            continue
        bt.display('tree', 'tree-'+str(line_number))

#bt.write_image(bt.GenerateHistogram())
bt.set_display_mode('matplotlib')
bt.display('tree')
#bt.write_image(bt.GenerateScatterplot())
#bt.write_image(bt.GenerateIncumbentPath())
#bt.write_image(bt.GenerateForecastImages())

bt.set_display_mode('file')
bt.display('tree', 'tree-final')
