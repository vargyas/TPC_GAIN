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



# set the extension of all figures
EXT = ["png","pdf","eps"]

# initial parameters (estimated based on past data)
init_mean = 400
init_width = 70
xmin = init_mean - 2.*init_width
xmax = init_mean + 2.*init_width
adc_pedestal = [2792, 2780, 2801, 2790]

from ROOT import TF1, TH1D, TH2F, TTree, TCanvas, TPad, gStyle, TFile, TLegend, TLine
import ROOT
import os
import warnings
from root_numpy import array2tree
from array import array
import numpy as np
import sys
import time

colors = [601, 634, 417, 433, 635, 619, 603, 636, 719, 435]

# Save time with not drawing the pop-up canvases
ROOT.gROOT.SetBatch(True)


def createFoilContour():
    """
    Draws original foil contour
    """

    xs1 = 78
    xs2 = 146
    xl1 = 56
    xl2 = 168
    l = []
    l.append( TLine(xs1, 2, xs2, 2) )
    l.append( TLine(xl1, 160, xl2, 160) )
    l.append( TLine(xl1, 160, xs1, 2) )
    l.append( TLine(xs2, 2, xl2, 160) )

    for il in l:
        il.SetLineColor(1)
        il.SetLineWidth(2)

    return l

def createSectorBound():
    x = [56, 112, 168]
    l=[]
    for ix in x:
        l.append( TLine( ix, 2, ix, 160) )
    for il in l:
        il.SetLineWidth(1)
        il.SetLineColor(1)
        il.SetLineStyle(9)

    return l

def createMapHisto(name):
    """
    Creates the mapping histogram
    """

    hMap = TH2F(name, "", 225, -0.5, 224.5, 161, -0.5, 160.5)
    return hMap

def processMap(inputFileName, outputFolderName):
    """
    Process ROOT TH2 histogram:
        - find mean and average value in each (x,y) point
        - plots the results
        - save fit parameters in the map to a tree
    """

    start = time.clock()

    # try to make directory
    try:
        os.makedirs(outputFolderName)
    except OSError:
        pass

    # Preparing to save data to a tree (using array2tree later)
    fit_array_tmp = list()

    # Read data from input file
    inFile = TFile(inputFileName, "READ")

    # Bind certain methods for speedup (python feature)
    hGainXY = inFile.Get("hGainXY")
    hGainDouble = inFile.Get("hDoubleCounted")
    #hGainXY.Add( hGainDouble )
    getBinContent = hGainXY.GetBinContent

    # Read gain peak histograms from file
    hGainZ = []
    igbin = 0
    while inFile.GetListOfKeys().Contains("hGainZ_{}".format(igbin)):
        hGainZ.append(inFile.Get("hGainZ_{}".format(igbin)))
        hGainZ[igbin].Sumw2()
        igbin += 1

    fitf = TF1("fitf", "gaus", 0, 1000)

    c = TCanvas("c","",800,600)

    for gbin in range(len(hGainZ)):
        # initialise, if fit succeeds will be assigned in the end of the loop
        #fit_array_tmp.append((ixb, iyb, 0.0, 0.0, 0.0, 0.0, 0.0, 0))
        #ic += 1

        ixb, iyb, izb = ROOT.Long(), ROOT.Long(), ROOT.Long()
        hGainXY.GetBinXYZ(gbin,ixb,iyb,izb)
        #gbin = getBin(ixb, iyb)

        # throw away small data
        peak_occ_init = getEntriesInRange( hGainZ[gbin], xmin+10, xmax-10 )
        #peak_occ_init = hGainZ[gbin].GetEntries()-hGainZ[gbin].GetBinContent(0)-hGainZ[gbin].GetBinContent(hGainZ[gbin].GetNbinsX())
        if peak_occ_init < 30:
            continue

        # throw away flat data
        #if hGainZ[gbin].GetStdDev() > 75.:
        #    continue

        # fit peak histogram
        fitf.SetParameters(0, init_mean, init_width)
        hGainZ[gbin].Fit(fitf, 'RQ')
        peak_mean = fitf.GetParameter(1)
        peak_width = fitf.GetParameter(2)
        chi2 = fitf.GetChisquare()
        ndf = fitf.GetNDF()

        #peak_mean = hGainZ[gbin].GetMean()
        #peak_width = hGainZ[gbin].GetStdDev()

        # throw away very suspicious fits
        if peak_mean < 0: continue
        if peak_mean > 2000: continue
        if peak_width <= 0: continue

        #hGainZ[gbin].Draw("")
        #leg = TLegend(0.7, 0.5, 0.9, 0.7)
        #leg.AddEntry(hGainZ[gbin], "[0],[1],[2]={:.3f}, {:.2f}, {:0.2f}".format(fitf.GetParameter(0),fitf.GetParameter(1),fitf.GetParameter(2)))
        #leg.AddEntry(hGainZ[gbin], "#chi^2/NDF={:.2f}/{:.0f}".format(chi2,ndf))
        #leg.Draw()
        #c.SaveAs('figs/{}/peakfinder/peak_X{}Y{}.{}'.format(folderName,ixb,iyb,EXT[0]))

        raw_occ  = getBinContent(gbin)
        peak_occ = getEntriesInRange(hGainZ[gbin], peak_mean-2.*peak_width, peak_mean+2.*peak_width)

        #fit_array_tmp[ic - 1] = (ixb, iyb, peak_mean, peak_width, raw_occ, peak_occ, chi2, ndf)
        fit_array_tmp.append( (ixb, iyb, peak_mean, peak_width, raw_occ, peak_occ, chi2, ndf) )


    stop = time.clock()
    print 'Generating map took [', stop - start, '] sec'

    fit_array = np.array(fit_array_tmp,
                         dtype=[('x', np.int),
                                ('y', np.int),
                                ('mean', np.float64),
                                ('width', np.float64),
                                ('raw_occ', np.float64),
                                ('peak_occ', np.float64),
                                ('chi2',np.float64),
                                ('ndf',np.int)
                                ])

    outFileName = outputFolderName+'/'+outputFolderName.split('/')[2]+'_tree.root'
    outFile = TFile(outFileName, "RECREATE")
    tree = array2tree(fit_array, "tree")
    outFile.Write()
    outFile.Close()

def drawMap(inputFolderName):
    """

    """

    infilename = inputFolderName+'/'+inputFolderName.split('/')[2]+'_tree.root'

    inFile = TFile.Open(infilename)
    dirName = os.path.dirname(infilename)
    tree = inFile.Get('tree')
    #tree.Print()

    # definition of a good cut
    goodCut = '(chi2/ndf<2)&&'
    goodCut+= '(width/mean<0.5)'

    c1 = TCanvas("c1", "c1", 1000, 1300)
    c1.Divide(2,3)
    gStyle.SetOptStat(0)
    #gStyle.SetOptTitle(0)

    c1.cd(1)
    tree.Draw("mean:x>>h1", goodCut+"&&mean>1", "scat")
    h1=tree.GetHistogram()
    h1.GetXaxis().SetRangeUser(0, 224)
    h1.GetYaxis().SetRangeUser(0, 1000)
    h1.SetMarkerColor(colors[1])
    h1.GetXaxis().SetTitle("x")
    h1.GetXaxis().SetTitle("mean")
    h1.Draw("scat")

    c1.cd(2)
    tree.Draw("mean:y>>h2", goodCut+"&&mean>1", "scat")
    h2=tree.GetHistogram()
    h2.GetXaxis().SetRangeUser(0, 160)
    h2.GetYaxis().SetRangeUser(0, 1000)
    h2.SetMarkerColor(colors[1])
    h2.GetXaxis().SetTitle("y")
    h2.GetXaxis().SetTitle("mean")
    h2.Draw("scat")


    c1.cd(3)
    tree.Draw("mean", goodCut+"&&mean>1", "")

    c1.cd(4)
    tree.Draw("width/mean",goodCut+"&&width/mean>0","")

    c1.cd(5)
    tree.Draw("raw_occ>>h3", goodCut+'&&(raw_occ>3)', "goff")
    h3=tree.GetHistogram()
    h3.GetXaxis().SetTitle("raw or peak occurence")
    h3.GetXaxis().SetRangeUser(3,1000)
    h3.SetLineColor(2)
    h3.Draw()
    tree.Draw("peak_occ>>h4",goodCut+'&&(peak_occ>3)',"goff")
    h4=tree.GetHistogram()
    h4.Draw("same")

    xmin = 10; xmax = 700
    efficiency = calcOccInRange(h4, xmin, xmax)/calcOccInRange(h3, xmin, xmax)*100.
    leg = TLegend(0.6,0.6,0.87,0.87)
    leg.SetBorderSize(0)
    leg.AddEntry(h3, "raw")
    leg.AddEntry(h4, "peak")
    leg.AddEntry(None, "{:.1f}%".format(efficiency), "")
    leg.Draw()
    c1.Update()

    c1.cd(6)
    tree.Draw("chi2/ndf","","")



    for iext in EXT: c1.SaveAs('{}/profiles.{}'.format(dirName,iext))

    tree.SetMarkerStyle(21)
    tree.SetMarkerSize(0.6)

    meanrange = '&&mean>150&&mean<550'
    tree.Draw("y:x:mean>>htemp",goodCut+"&&mean>1","colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "mean")

    tree.Draw("y:x:(width/mean)>>htemp", goodCut+"&&(width/mean)>0", "colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "width")

    tree.Draw("y:x:raw_occ>>htemp", "raw_occ>1", "colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "raw_occ")

    tree.Draw("y:x:peak_occ>>htemp", goodCut+'&&(peak_occ>1)', "colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "peak_occ")

    tree.Draw("y:x:chi2/ndf>>htemp", goodCut, "colz")
    h = tree.GetHistogram()
    drawTreeMap(h, dirName, "chi2ndf")


def getEntriesInRange(h, xmin, xmax):
    xbinmin = h.FindBin(xmin)
    xbinmax = h.FindBin(xmax)
    return h.Integral(xbinmin, xbinmax)


def drawTreeMap(h, dir, name):
    """
    Draws and saves given histogram
    :param h: histogram to save
    :param name: save name
    :return:
    """
    c1 = TCanvas("c1", "c1", 1000, 1300)
    ROOT.gPad.SetRightMargin(0.15)
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)

    h.GetXaxis().SetRangeUser(50., 172.)
    h.GetYaxis().SetRangeUser(0., 161.)
    h.Draw("colz")
    l = createFoilContour()
    for il in l: il.Draw()
    ls = createSectorBound()
    for il in ls: il.Draw()
    for iext in EXT: c1.SaveAs('{}/{}.{}'.format(dir,name,iext))


def calcOccInRange(h_occ, xmin, xmax):
    ixmin = h_occ.FindBin(xmin)
    ixmax = h_occ.FindBin(xmax)
    h_occ_GetX = h_occ.GetBinCenter
    h_occ_GetY = h_occ.GetBinContent
    all = 0
    for ib in range(ixmin, ixmax):
        all += h_occ_GetX(ib) * h_occ_GetY(ib)
    return all


### M A I N   P R O G R A M
#input is foldername + run name
#example IROC_14_a 127

inputFileName = os.path.normpath(os.path.join(sys.argv[1],'data',sys.argv[2]))
outputFolderName = os.path.normpath(os.path.join(sys.argv[1],'figs',sys.argv[2])).replace('GemQa2','')
outputFolderName = outputFolderName.replace('.root','')

print inputFileName
print outputFolderName


print ""
#processMap(inputFileName, outputFolderName)
print ""
drawMap(outputFolderName)
print ""
