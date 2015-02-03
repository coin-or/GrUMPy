from src.grumpy import BBTree
import sys
import StringIO

try:
    from PIL import Image as PIL_Image
except ImportError:
    PIL_INSTALLED = False
else:
    PIL_INSTALLED = True

if __name__ == '__main__':
    bt = BBTree()
    bt.set_display_mode('pygame')
    line_number = 0
#    file_ = open('air04.bak', 'r')
    for line in sys.stdin:
#    for line in file_:
        bt.ProcessLine(line)
        print 'line', line_number, 'processed.'
        line_number = line_number+1
        if line_number%10000 != 0:
            continue
#        imagefile = open('tree-'+str(line_number)+'.png','w')
#        imagefile.write(bt.GenerateTreeImage())
#        imagefile.close()
        if bt.root is not None:
#            gnuplot_image = bt.GenerateHistogram()
            gnuplot_image = StringIO.StringIO(bt.GenerateTreeImage())
#            gnuplot_image = bt.GenerateScatterplot()
#            gnuplot_image = bt.GenerateIncumbentPath()
#            gnuplot_image = bt.GenerateForecastImages()
            if gnuplot_image is not None and PIL_INSTALLED:
                im = PIL_Image.open(gnuplot_image)
                im.show()
#            bt.display_all()

    gnuplot_image = StringIO.StringIO(bt.GenerateTreeImage())
    if gnuplot_image is not None and PIL_INSTALLED:
        im = PIL_Image.open(gnuplot_image)
        im.show()
#   bt.display_all()
#    raw_input("Press a key to exit")
