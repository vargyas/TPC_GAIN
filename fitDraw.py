###################################################################
### This macro reads the gain scan output and performs two tasks.
### 1. Saves the preprocessed output to a ROOT tree.
### 2. Plots spectrum, and map
###
### Author: Vargyas, Marton
### Email: mvargyas@cern.ch
### 2016
###################################################################

# to run on all data (depreciated):
# for f in data/*root; do tmpName=${f##*/}; python fitDraw.py ${tmpName/.root/} ; done

# set the extensions of all figures
EXT = ["png","pdf","eps"]

# initial parameters (estimated based on past data)
init_mean = 400
init_width = 70
xmin = init_mean - 2*init_width
xmax = init_mean + 2*init_width
adc_pedestal = [2792.0, 2780.0, 2801.0, 2790.0]

import os
import numpy as np
import sys
import time
from scipy.stats import norm
from scipy.stats import chisquare
import matplotlib.pyplot as plt
import matplotlib.cm as cm


colors = [601, 634, 417, 433, 635, 619, 603, 636, 719, 435]

def processMap(inputFileName):
    """
    - find mean and average value in each (x,y) point by Gaussian fit
    - plots the results
    - save fit parameters in the map to a tree
    """

    # control whether to plot Gaussian fits
    plotFitProcess = True

    start = time.clock()
    #outputFolderName = os.path.dirname(os.path.realpath(inputFileName))

    fit_array_tmp = list()

    # Read data from input file
    inFile = open(inputFileName, "r")
    lines = inFile.readlines()
    nl = len(lines)

    start = time.clock()
    ic = 0
    for il in range(nl):
        words = lines[il].split()
        x = int(words[0])
        y = int(words[1])
        gbin = int(words[2])
        adc = words[3:]
        mu = -1
        std = -1
        nraw = len(adc)
        # Fill with garbage to have all points
        fit_array_tmp.append((x, y, mu, std, nraw, 0.0, 0.0, 0))

        if len(adc)>30:
            adcf = map(float, adc)
            mu = np.mean(adcf)
            std= np.std(adcf)
            # ESTIMATE RESULTS AND SLICE THE ARRAY
            # OTHERWISE THE FIT PICKS UP ON OUTLIERS
            xmin = mu - 2 * std
            xmax = mu + 2 * std

            adcf_sliced = [_ix for _ix in adcf if _ix>xmin and _ix<xmax]
            # Fit a normal distribution to the data:
            mu, std = norm.fit(adcf_sliced)
        
            ndf = 1
            xmin = mu-2*std
            xmax = mu+2*std
            # find occurence of values that are in +-2sigma off of the mean
            npeak = 0
            for iadc in adcf_sliced:
                if iadc>xmin and iadc<xmax:
                    npeak+=1

            fit_array_tmp[ic] = (x, y, mu, std, nraw, npeak, 0, ndf)
            
            if plotFitProcess:
                plt.clf()

                count, bins, ignored = plt.hist(adcf_sliced, 30, normed=True)
                plt.plot(bins, 1/(std * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * std**2) ), linewidth=2, color='r')
                #p = norm.pdf(plot_x, mu, std)
                #plt.plot(plot_x, p, 'k', linewidth=2)
                title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
                plt.title(title)
                plt.savefig("fig_fit/fit_X{}_Y{}.png".format(x, y))
                #plt.show()

        ic += 1

    pass
    stop = time.clock()
    print 'Generating map took [', stop - start, '] sec'

    fit_array = np.array(fit_array_tmp,
                         dtype=[('x', np.int),
                                ('y', np.int),
                                ('mean', np.float64),
                                ('sigma', np.float64),
                                ('nraw', np.float64),
                                ('npeak', np.float64),
                                ('chi2',np.float64),
                                ('ndf',np.float64)
                                ])

    outputFileName = inputFileName.replace('PREMAP.txt','MAP.txt')
    np.savetxt(outputFileName, fit_array)

    return outputFileName



# DRAW MAP:

def drawMap(inputFileName):

    X = np.linspace(0,224,224)
    Y = np.linspace(0,160,160)
    XY = np.linspace(0, 224*160, 224*160)


    array = np.loadtxt(inputFileName, dtype={'names': ('x', 'y', 'mean', 'sigma', 'nraw', 'npeak', 'chi2', 'ndf'),
                                   'formats': (float, float, float, float, float, float, float, float)})



    grid = array['npeak'].reshape(224, 160).T

    #plt.figure(figsize=(100,100))
    #plt.imshow(grid, cmap=cm.rainbow, vmin=30,vmax=1600)
    #plt.clim(50, 1600)
    #plt.gca().invert_yaxis()
    #plt.colorbar()
    #plt.show()

    hist=np.histogram(array['mean'])
    plt.hist(array['mean'], bins=np.arange(0,550))
    plt.show()

    '''
    plt.scatter(array['x'], array['mean'], c = array['chi2'], cmap = cm.hsv)
    plt.ylim(0, 1000)
    plt.show()

    ### create profile multifig
    grid_size = (4, 2)

    plt.subplot2grid(grid_size, (0,0), rowspan=3,colspan=1)
    #plt.ylim(0, 1000)
    plt.plot(array['nraw'])

    plt.subplot2grid(grid_size, (0, 1), rowspan=3, colspan=1)
    # plt.ylim(0, 1000)
    plt.plot(array['npeak'])

    plt.tight_layout()
    plt.show()
    '''
    #plt.contour(array['x'], array['y'], array['mean'])
    #plt.show()
"""
    tree = array2tree(array, 'tree')

    tree.Print()
    dirName = os.path.dirname(inputFileName)
    # definition of a good cut
    goodCut = ''



    c1 = TCanvas("c1", "c1", 1000, 1300)
    c1.Divide(2,3)
    gStyle.SetOptStat(0)
    #gStyle.SetOptTitle(0)

    c1.cd(1)
    tree.Draw("mean:x>>h1", "mean>1", "scat")
    h1=tree.GetHistogram()
    h1.GetXaxis().SetRangeUser(0, 224)
    h1.GetYaxis().SetRangeUser(0, 1000)
    h1.SetMarkerColor(colors[1])
    h1.GetXaxis().SetTitle("x")
    h1.GetXaxis().SetTitle("mean")
    h1.Draw("scat")

    c1.cd(2)
    tree.Draw("mean:y>>h2", "mean>1", "scat")
    h2=tree.GetHistogram()
    h2.GetXaxis().SetRangeUser(0, 160)
    h2.GetYaxis().SetRangeUser(0, 1000)
    h2.SetMarkerColor(colors[1])
    h2.GetXaxis().SetTitle("y")
    h2.GetXaxis().SetTitle("mean")
    h2.Draw("scat")

    c1.cd(3)
    tree.Draw("mean", "", "")

    c1.cd(4)
    tree.Draw("sigma/mean","","")

    c1.cd(5)
    tree.Draw("nraw>>h3", "nraw>0", "goff")
    h3=tree.GetHistogram()
    h3.GetXaxis().SetTitle("raw or peak occurence")
    h3.GetXaxis().SetRangeUser(3,1000)
    h3.SetLineColor(2)
    h3.Draw()
    tree.Draw("npeak>>h4","npeak>0","goff")
    h4=tree.GetHistogram()
    h4.Draw("same")

    leg = TLegend(0.6,0.6,0.87,0.87)
    leg.SetBorderSize(0)
    leg.AddEntry(h3, "raw")
    leg.AddEntry(h4, "peak")
    leg.Draw()
    c1.Update()

    c1.cd(6)
    tree.Draw("chi2/ndf","","")



    for iext in EXT: c1.SaveAs('{}/profiles.{}'.format(dirName,iext))

    tree.SetMarkerStyle(21)
    tree.SetMarkerSize(0.6)

    tree.Draw("y:x:mean>>htemp","","colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "mean")

    tree.Draw("y:x:(sigma/mean)>>htemp", "sigma/mean>0", "colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "width")

    tree.Draw("y:x:nraw>>htemp", "", "colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "raw_occ")

    tree.Draw("y:x:npeak>>htemp", "", "colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "peak_occ")

    tree.Draw("y:x:chi2/ndf>>htemp", "chi2>0", "colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "chi2ndf")

"""

"""
def drawTreeMap(h, dir, name):

    Draws and saves given histogram
    :param h: histogram to save
    :param name: save name
    :return:

    c1 = TCanvas("c1", "c1", 1000, 1300)
    ROOT.gPad.SetRightMargin(0.15)
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)

    h.GetXaxis().SetRangeUser(50., 172.)
    h.GetYaxis().SetRangeUser(0., 161.)
    h.Draw("colz")
    #l = createFoilContour()
    #for il in l: il.Draw()
    #ls = createSectorBound()
    #for il in ls: il.Draw()
    for iext in EXT: c1.SaveAs('{}/{}.{}'.format(dir,name,iext))

"""

### M A I N   P R O G R A M
# input is mapfile.root
# example: python fitDraw.py IROC_14_a/Run127/Run127PREMAP.txt


print ""
outname = processMap(sys.argv[1])
print ""
drawMap(outname)
#print ""
