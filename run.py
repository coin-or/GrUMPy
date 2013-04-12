from baktree import BAKTree
import sys

#
# tree
# histogramuare
# scatterplot
# imcumbentpath
#

#
# GenerateTreeImage done
# GenerateHistogram done
# GenerateScatterplot (we have very few nodes to draw) done
# GenerateIncumbentPath
# GenerateAllIncumbentPaths
# CreateAnimatedImages (not relevant)
#
# GeneratePredictionImages
# measures.png irrelevant, drawn at the end of the run.
#

if __name__ == '__main__':
    bt = BAKTree()
    bt.set_display_mode('pygame')
    line_number = 0
    file_ = open('p0201_GLPK.in', 'r')
#    for line in sys.stdin:
    for line in file_:
        bt.ProcessLine(line)
        print 'line', line_number, 'processed.'
        line_number = line_number+1
        if line_number%200 != 0:
            continue
        bt.write_image(bt.GenerateTreeImage(), filename = "image"+str(line_number)+".png")
#        if bt.root is not None:
#            gnuplot_image = bt.GenerateHistogram()
#            gnuplot_image = bt.GenerateTreeImage()
#            gnuplot_image = bt.GenerateScatterplot()
#            gnuplot_image = bt.GenerateIncumbentPath()
#            gnuplot_image = bt.GenerateForecastImages()
#            if gnuplot_image is not None:
#                bt.display_image(gnuplot_image)
#            bt.display_all()
    
    bt.write_image(bt.GenerateTreeImage())

