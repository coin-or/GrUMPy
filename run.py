from baktree import BAKTree
import sys

#
# tree
# histogram
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
#    file_ = open('p0201_GLPK.in', 'r')
    for line in sys.stdin:
#    for line in file_:
        bt.ProcessLine(line)
        print 'line', line_number, 'processed.'
        line_number = line_number+1
        if line_number%100 != 0:
            continue
        if bt.root is not None:
#            gnuplot_image = bt.GenerateHistogram()
            gnuplot_image = bt.GenerateTreeImage()
#            gnuplot_image = bt.GenerateScatterplot()
#            gnuplot_image = bt.GeneratePredictionImages()
            if gnuplot_image is not None:
                bt.display_image(gnuplot_image)
