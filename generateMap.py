###################################################################
### This macro reads the gain scan output (.ebe) and converts it to
### an X-Y map ROOT file in the same directory, with the same name
### as the input, along with the adc values in each point.
###
### Author: Vargyas, Marton
### Email: mvargyas@cern.ch
### 2016
###################################################################


# initial parameters (estimated based on past data)
init_mean = 440
init_width = 100
adc_pedestal = [2792, 2780, 2801, 2790]

# from ROOT import TF1, TH1D, TH2F, TFile
#import ROOT
from array import array
import numpy as np
import sys
import time


# save a little time with not drawing the canvas
# ROOT.gROOT.SetBatch(True)


def getGlobalIndex(x, y):
    return 224*y + x


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

    aADC = []
    for i in range(224 * 160):
        aADC.append(array('f'))

    with open(fileName) as f:
        lines = f.readlines()
        nl = len(lines)

        start = time.clock()
        for il in range(nl):

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
            # if x ==-2 and y ==-2:
            #    hDoubleCountedFill(xl, yr)
            #    hDoubleCountedFill(xr+112, yr)

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

            x = int(x)
            y = int(y)
            gbin = getGlobalIndex(x, y)
            aADC[gbin].append(adc)

        stop = time.clock()
        print '\nProcessing file took [', stop - start, '] sec'

        outFileName = fileName.replace('.ebe', 'PREMAP.txt')
        outFileName = outFileName.replace('GemQa2','')
        print 'Saving map to: {}'.format(outFileName)

        outF = open(outFileName, 'w+')
        for ix in range(224):
            for iy in range(160):

                gbin = getGlobalIndex(ix, iy)
                adc_str = ""
                if len(aADC[gbin])>30:
                    for iadc in aADC[gbin]:
                        adc_str+=str(iadc)+"\t"

                outF.writelines( ["{:d}\t{:d}\t{:d}\t{:s}\n".format(ix, iy, getGlobalIndex(ix,iy), adc_str)] )

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
    Find cluster (=continuous array)
    from short array
    """

    mean = arr[0]
    # break if it is not continuous
    for i in range(len(arr) - 1):
        mean += arr[i + 1]
        if arr[i]+1 != arr[i + 1]:
            return -1
    mean = float(mean) / float(len(arr))
    return mean


def findCluster(xl, xr, yl, yr):
    """
    Find the valid cluster x or y
    > 0: valid, left of right
    -1: not valid, only one valid hit
    -2: not valid, two valid hits
    """

    x = -1
    y = -1

    # LEFT
    if xl > 0 and xr < 0:
        if yl > 0 and yr < 0:
            x = xl
            y = yl
    # RIGHT
    elif xl < 0 and xr > 0:
        if yl < 0 and yr > 0:
            x = xr + 112
            y = yr
    elif xl > 0 and xr > 0:
        x = -2
        if yl > 0 and yr > 0:
            y = -2
        else:
            y = -1
    else:
        x = -1
        if yl > 0 and yr > 0:
            y = -2
        else:
            y = -1

    return x, y


### M A I N   P R O G R A M
print ""
generateMap(sys.argv[1])
print ""
