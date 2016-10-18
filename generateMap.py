###################################################################
### This macro reads the gain scan output and converts it to
### an X-Y map, along with the adc values in each point.
###
### Input is the filename with path, creates the map in that folder
###
### Author: Vargyas, Marton
### Email: mvargyas@cern.ch
### 2016
###################################################################


# initial parameters (estimated based on past data)
init_mean = 440
init_width = 100
adc_pedestal = [2792, 2780, 2801, 2790]

from ROOT import TF1, TH1D, TH2F, TFile
import ROOT
from array import array
import sys
import time

# save a little time with not drawing the canvas
ROOT.gROOT.SetBatch(True)



def createMapHisto(name):
    """
    Creates the mapping histogram
    """

    hMap = TH2F(name, "", 225, -0.5, 224.5, 161, -0.5, 160.5)
    # hMap = TH2F(name,"",112, 55.5, 167.5,160, -0.5, 159.5)
    return hMap


def getNextIndex(words, prev_index):
    """
    Reads x or y index and gets the next index
    (order is xl, yl, xr, yr)
    """

    index_len = int(words[prev_index])
    if index_len == 0:
        next_index = prev_index + 1
    else:
        next_index = prev_index + index_len + 1
    return next_index


def generateMap(fileName):
    """
    Load text file sent by the raspberry pi
    into a ROOT TH2 histogram with basic pre-processing:
        - rejecting separate clusters in one event
        - save results to ROOT file
    """

    print 'generating map from: ', fileName
    outFileName = fileName.replace('.ebe', '.root')
    outFile = TFile(outFileName, "RECREATE")

    hGainXY = createMapHisto("hGainXY")
    hDoubleCounted = createMapHisto("hDoubleCounted")
    # hGainXY.Rebin2D(2,2)
    # hGainXY.RebinX(2); hGainXY.RebinY(2)

    aGainADC = []
    hGainZ = []

    for igbin in range((hGainXY.GetNbinsX()+1) * (hGainXY.GetNbinsY()+1)):
        aGainADC.append(array('f'))
        hGainZ.append(TH1D("hGainZ_{}".format(igbin), "", 1000, 0.5, 1000.5))

    getXbin = hGainXY.GetXaxis().FindBin
    getYbin = hGainXY.GetYaxis().FindBin
    hGainXYFill = hGainXY.Fill
    hDoubleCountedFill = hDoubleCounted.Fill
    getBin = hGainXY.GetBin

    with open(fileName) as f:
        lines = f.readlines()
        print 'Processing file: {}'.format(fileName)
        nl = len(lines)

        # printProgress(0, nl, prefix = loadStr, suffix = 'complete', barLength = 50)

        start = time.clock()
        for il in range(nl):

            # printProgress(il, nl, prefix = loadStr, suffix = 'complete', barLength = 50)

            xarrl = []
            yarrl = []
            xarrr = []
            yarrr = []
            words = lines[il].split()
            ixl_index = 3
            iyl_index = getNextIndex(words, ixl_index)
            ixr_index = getNextIndex(words, iyl_index)
            iyr_index = getNextIndex(words, ixr_index)
            iyr_index_up = getNextIndex(words, iyr_index)

            # collect x left
            for ixl in range(ixl_index + 1, iyl_index):
                xarrl.append(int(words[ixl]))
            # collect y left
            for iyl in range(iyl_index + 1, ixr_index):
                yarrl.append(int(words[iyl]))
            # collect x right
            for ixr in range(ixr_index + 1, iyr_index):
                xarrr.append(int(words[ixr]))
            # collect y right
            for iyr in range(iyr_index + 1, iyr_index_up):
                yarrr.append(int(words[iyr]))

            # it rejects empty channel arrays to save time, then
            # it rejects non-continuous channels
            if len(xarrl) > 0:
                xl = findClusterPart(xarrl)
            else:
                xl = -1
            if len(xarrr) > 0:
                xr = findClusterPart(xarrr)
            else:
                xr = -1
            if len(yarrl) > 0:
                yl = findClusterPart(yarrl)
            else:
                yl = -1
            if len(yarrr) > 0:
                yr = findClusterPart(yarrr)
            else:
                yr = -1

            # we accept only one cluster pair (x-y)
            # decide here which will be kept
            x, y = findCluster(xl, xr, yl, yr)

            # fill the double counted histo
            if x ==-2 and y ==-2:
                hDoubleCountedFill(xl, yr)
                hDoubleCountedFill(xr+112, yr)

            # only fill if both x and y are valid clusters
            if x < 0 or y < 0:
                continue

            # get ADC index and value
            adc_all = words[iyr_index_up + 1:iyr_index_up + 1 + 4]
            which_adc = findAdcIndex(x)
            adc = float(adc_all[which_adc])
            # for some reason half value is stored, correcting here
            if which_adc == 1:
                adc = 2. * adc

            # subtract pedestal after the correction
            adc = adc - adc_pedestal[which_adc]

            # also this correction in the awk file
            if which_adc == 0:
                adc = 1.4 * adc

            xbin = getXbin(x)
            ybin = getYbin(y)
            hGainXYFill(x, y)

            gbin = getBin(xbin, ybin) # gbin starts from 0!
            hGainZ[gbin].Fill(adc)
            aGainADC[gbin].append(adc)
            pass

    outFile.Write()
    outFile.Close()

    stop = time.clock()
    print '\nProcessing file took [', stop - start, '] sec'


def findAdcIndex(x):
    """
    Find the Adc index from x value
    """

    if x >= 0 and x < 56:
        return 0
    elif x >= 56 and x < 112:
        return 1
    elif x >= 112 and x < 168:
        return 2
    elif x >= 168 and x < 224:
        return 3


def findClusterPart(arr):
    """
    Find cluster from short array
    """

    mean = arr[0]
    # break if it is not continuous
    for i in range(len(arr) - 1):
        mean += arr[i + 1]
        if arr[i] != arr[i + 1] - 1:
            return -1
    mean = float(mean) / float(len(arr))
    return mean


def findCluster(xl, xr, yl, yr):
    """
    Find the valid cluster x or y
    -1: valid, only one valid hit
    -2: not valid, two valid hits
    """

    if xl > 0 and xr < 0:
        x = xl
    elif xl < 0 and xr > 0:
        x = xr + 112
    elif xl > 0 and xr > 0:
        x = -2
    else:
        x = -1

    if yl > 0 and yr < 0:
        y = yl
    elif yl < 0 and yr > 0:
        y = yr
    elif yl > 0 and yr > 0:
        y = -2
    else:
        y = -1

    return x, y





### M A I N   P R O G R A M
print ""
generateMap(sys.argv[1])
print ""
