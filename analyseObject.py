#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

from re import search

from gPhoton.gAperture import gAperture

from numpy import mean, median, argmax, genfromtxt, concatenate, where
from math import sqrt, floor, ceil, log10

import numpy as np
import numpy.ma as ma

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.font_manager as font_manager
import matplotlib

from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

from interactivePlotting import AnnoteFinder

from bisect import bisect_right,bisect_left

from astropy.constants import L_sun, R_sun, sigma_sb
from astropy import units as u
from astropy.table import Table
from astropy.io import fits, ascii
from astropy.time import Time

import os
import subprocess

from Queue import Queue
from threading import Thread
from time import time

from textwrap import wrap

# Global constants

Mbol_sun = 4.74 # Ryden cosmology book

galexEqWth = 790 # angstroms, crude est from Rachel using weird table and IDL
galexBPpercentE = 0.096 # using Osten 2015 and IDL plank integration (redo)

dataPath = '/Users/cbrasseur/Documents/FlareWork/data/'
galexLCPath = dataPath + 'flareFiles/galexLcs/'
keplerLCPath = dataPath + 'flareFiles/keplerLcs/'

fontpath = '/Library/Fonts/Georgia.ttf'

def initLogfile(LFPNT):
    LFPNT.write("GALEX Object ID\n")
    LFPNT.write("ra,dec\n")
    LFPNT.write("Light Curve csv file(s)\n")
    LFPNT.write("Anything else\n")


def getObsIntervals(intLst):

    intervals = []

    for intv in intLst:
        intSearch = search('\((.*):(.*)\)',intv)
        intervals.append([float(intSearch.group(1)),
                          float(intSearch.group(2))])       
    return intervals

                         

def makeLCFiles(kgObjStr): 

    objLst = kgObjStr.strip('\n').split(',')

    # open the logfile (APPEND mode)
    LOG = open('KGLightcurveLog.txt','a')

    if not LOG.tell():
        initLogfile(LOG)

    LOG.write('\n' + objLst[1] + '\n')
    LOG.write(objLst[4] + ',' + objLst[5] + '\n')


    # Make light curves for each time interval
    

    timeIntervals = getObsIntervals(objLst[6:])

    for obsInt in timeIntervals:
        csvFle = objLst[1]+"_" + str(obsInt[0]) + "_" + str(obsInt[1]) + ".csv"

        ltCv = gAperture(band='NUV',
                         skypos=[float(objLst[4]),float(objLst[5])],
                         trange=obsInt,
                         radius=0.004,
                         annulus=[.005,.007],
                         stepsz=10,
                         csvfile=csvFle)
        #print("done one")

        LOG.write(csvFle + '\n')
        

def getCombinedLC(kgObjStr):

    objLst = kgObjStr.strip('\n').split(',')

    timeIntervals = getObsIntervals(objLst[6:])

    ltCv = gAperture(band='NUV',
                     skypos=[float(objLst[4]),float(objLst[5])],
                     tranges=timeIntervals,
                     radius=0.004,
                     annulus=[.005,.007],
                     stepsz=10)
    return ltCv


# statistic for flare:
# max value is 3.5 sigma above mean
# at least one point directly to right or left
# is 2 sigma above mean
def isFlare(valLst):
    maxLoc = argmax(valLst)

    A = np.mean(valLst)
    
    M = valLst[maxLoc] #maxVal
    sigM = sqrt(M)
    
    if ((M-A)/sigM) < 3.5: # maxVal is less than 3 sigma from mean
        return False

    # look at value after maxVal (unless maxVal is last value)
    if maxLoc != (len(valLst) - 1):
        M = valLst[maxLoc+1]
        sigM  = sqrt(M)

        if ((M-A)/sigM) >= 2: # value next to maxVal is 2 sigma or above
            return True
        
    # look at value before maxVal (unless maxVal is first value)
    if maxLoc != 0:       
        M = valLst[maxLoc-1]
        sigM = sqrt(M)

        if ((M-A)/sigM) >= 2: # value next to maxVal is 2 sigma or above
            return True

    # values to both sides of maxVal are < 2 sigma
    return False



def findFlareRough(kgObjStr, quiet = True):

    objLst = kgObjStr.strip('\n').split(',')

    timeIntervals = getObsIntervals(objLst[6:])

    lcfilename = objLst[1] + "_LC.csv"
    #print('gAperture('+'NUV'+','+'[' + objLst[4] +','+ objLst[5] + ']' + ',' + str(timeIntervals) + ',0.004,[.005,.007],10,'+lcfilename + ')')
    #return
    ltCv = gAperture(band='NUV',
                     skypos=[float(objLst[4]),float(objLst[5])],
                     tranges=timeIntervals,
                     radius=0.004,
                     annulus=[.005,.007],
                     stepsz=10,
                     csvfile = lcfilename)

    
    
    cnts = ltCv['counts']
    
    if isFlare(cnts):  # FLARE
        if not quiet:
            print("GALEX object", objLst[1], "FLARES. Filename: ", lcfilename)

    else: # no flare

        if not quiet:
            print("GALEX object", objLst[1], "does not flare.")
        os.remove(lcfilename)
        lcfilename = ""

    # log results
    with open('KGFlareResults.csv','a') as LOG:

        if not LOG.tell():
            LOG.write("galex_id,ra,dec,min_cnt,max_cnt,ave_cnt,lc_filename\n")

        lineLst = [objLst[1],objLst[4],objLst[5],str(min(cnts)),str(max(cnts)),str(np.mean(cnts)),lcfilename]

        LOG.write(','.join(lineLst) + '\n')
    
    return lcfilename

   
    

def firstPassFindFlares(startLine = 0, lines = 1, quiet = True):

    ts = time()
    
    # not most efficient, but putting the whole file in an array
    with open("KGSimultaneousObservations.csv",'r') as KGOBJ:
        objStrings = KGOBJ.readlines()[1:] # first line is col headers

    if (lines == 'all') or ((startLine + lines) > len(objStrings)):
        endLine = len(objStrings)
    else:
        endLine = startLine + lines

    flareLCs = []

    while (startLine < endLine):
        lc = findFlareRough(objStrings[startLine],quiet)
        if not lc == "":
            flareLCs.append(lc)

        startLine +=1

    print("Elapsed time in seconds:",time()-ts)

    return flareLCs

class FindFlareWorker(Thread):
    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            kgObjStr, quiet = self.queue.get()
            findFlareRough(kgObjStr, quiet)
            self.queue.task_done()

            
def firstPassFindFlares_withThreading(startLine = 0, lines = 1, quiet = True):

    ts = time()
    # not most efficient, but putting the whole file in an array
    with open("KGSimultaneousObservations.csv",'r') as KGOBJ:
        objStrings = KGOBJ.readlines()[1:] # first line is col headers

    if (lines == 'all') or ((startLine + lines) > len(objStrings)):
        endLine = len(objStrings)
    else:
        endLine = startLine + lines

    flareLCs = []

    # Create a queue to communicate with the worker threads
    queue = Queue()

    # Create 6 worker threads
    for x in range(6):
        worker = FindFlareWorker(queue)
        worker.start()

    # Put the tasks into the queue as a tuple
    for i in range(startLine,endLine):
        queue.put((objStrings[i], False))

    queue.join()
    
    print("Elapsed time in seconds:",time()-ts)

    

def findIntervals(galexID):

    grepRes = subprocess.Popen(['grep',galexID,dataPath + "KGSimultaneousObservations.csv"],stdout=subprocess.PIPE)
    obsStr = grepRes.stdout.read().split('\n')[0] # if duplicates only want first copy
    
    if obsStr == '':
        print("No object string...")
        return []

    obsLst = obsStr.split(',')

    return getObsIntervals(obsLst[6:])
    

    
        
def graphLightCurves(LCfilenames):
    if not type(LCfilenames) == list:
        LCfilenames = [LCfilenames]

    for LCdatafile in LCfilenames:

        # get full light curve
        data = genfromtxt(LCdatafile,delimiter=',',names=True)

        # get galex time intervals
        obsInts = findIntervals(LCdatafile[-26:-7])

        # make the tables for graphing
        countByTime=[]
        ii = 0
        for start,stop in obsInts:
            fi = bisect_right(data['t_mean'],stop)
            countByTime.append((data['t_mean'][ii:fi],data['counts'][ii:fi]))
            ii = fi

        # set up the subplots
        nplots = len(countByTime)
        if not nplots:
            print("Find intervals says there are not plots to make.")
            continue
        ncols = floor(sqrt(nplots))
        nrows = ceil(nplots/ncols)

        # for sharing y axis limits
        ax1 = None

        for i in range(nplots):
            if not ax1:
                ax1 = plt.subplot(nrows, ncols, i+1)
            else:
                plt.subplot(nrows, ncols, i+1, sharey=ax1)

            error = [sqrt(x) for x in countByTime[i][1]]
                        
            plt.errorbar(countByTime[i][0], countByTime[i][1], yerr=error,fmt='b-o')
            plt.axhline(median(countByTime[i][1]), color='g')
            plt.xlabel('t_mean')
            plt.ylabel('counts (poisson error)')

        plt.suptitle('GALEX object id: ' + LCdatafile[:-7])
        
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())

        plt.show()


# statistic for flare:
# max value is 3.5 sigma above median
# at least one point directly to right or left
# is 2 sigma above median
# max value is further from median than min value
def isReallyFlare(valLst, med = None):
    maxLoc = argmax(valLst)
    minVal = min(valLst)

    A = median(valLst)
    
    M = valLst[maxLoc] #maxVal
    sigM = sqrt(M)

    if (M - A) < (A - minVal): # maxVal is further from median than minVal
        return False
    
    if ((M-A)/sigM) < 3.5: # maxVal is less than 3 sigma from mean
        return False

    # look at value after maxVal (unless maxVal is last value)
    if maxLoc != (len(valLst) - 1):
        M = valLst[maxLoc+1]
        sigM  = sqrt(M)

        if ((M-A)/sigM) >= 2: # value next to maxVal is 2 sigma or above
            return True
        
    # look at value afer maxVal (unless maxVal is first value)
    if maxLoc != 0:       
        M = valLst[maxLoc-1]
        sigM = sqrt(M)

        if ((M-A)/sigM) >= 2: # value next to maxVal is 2 sigma or above
            return True

    # values to both sides of maxVal are < 2 sigma
    return False


    
# applied second round of tests on files that passed the rough find test
def secondPassTest(LCfilenames, quiet = True):
    if not type(LCfilenames) == list:
        LCfilenames = [LCfilenames]
        
    for LCdatafile in LCfilenames:

        # get full light curve
        data = genfromtxt(LCdatafile,delimiter=',',names=True)

        
        if isReallyFlare(data['counts']):  # FLARE
            if not quiet:
                print("GALEX object", LCdatafile[:-7], "FLARES.")
            os.rename(LCdatafile,"./flareFiles/galexLcs/"+LCdatafile)

        else: # no flare

            if not quiet:
                print("GALEX object",LCdatafile[:-7], "does not flare.")
            os.rename(LCdatafile,"./noFlareFiles/"+LCdatafile)

            
def findGaps(endl = None):

    # getting the galex IDs from the simultaneous list
    with open("KGSimultaneousObservations.csv",'r') as KGOBJ:
        objStrings = KGOBJ.readlines()[1:] # first line is col headers
    gIDs2Obs = [x.strip('\n').split(',')[1] for x in objStrings]

    # getting the galex IDs from the analysed list
    with open("KGFlareResults.csv",'r') as KGANA:
        objStrings = KGANA.readlines()[1:] # first line is col headers
    gIDsAnal = [x.strip('\n').split(',')[0] for x in objStrings]

    lineNosToObs = []

    if (not endl) or (endl > len(gIDs2Obs)):
        endl = len(gIDs2Obs)
        

def findFlareEdges(fluxLst, flarePeakInd, fluxMed):
    flareBeg = flarePeakInd
    flareEnd = flarePeakInd
    
    while (fluxLst[flareBeg] > fluxMed) and (flareBeg > 0):
        flareBeg -= 1
        
    while (fluxLst[flareEnd] > fluxMed) and (flareEnd < (len(fluxLst)-1)):
        flareEnd += 1

    return flareBeg,flareEnd
        

def getFlareWidths(times,fluxes,counts,maxGap):
    
    fluxMed = median(fluxes)
    
    errors = [sqrt(x)*(y/x) if not x == 0 else 0 for x,y in zip(counts,fluxes)]
    

    # we only consider the light curve continuous when the gap between
    # flux measurements is <= maxGap (in seconds) so we want to divide
    # up the light curve into continuous segments
    timesFluxesErrors = []
    medAndErrs = []
    startInt = 0
    for endInt in range(1,len(times)):
        if (times[endInt] - times[endInt-1]) > maxGap:
            timesFluxesErrors.append((times[startInt:endInt],
                                      fluxes[startInt:endInt],
                                      errors[startInt:endInt]))
            fm = median(fluxes[startInt:endInt])
            cm = median(counts[startInt:endInt])
            medAndErrs.append((fm,sqrt(cm)*(fm/cm)))
            startInt = endInt
    timesFluxesErrors.append((times[startInt:],
                                 fluxes[startInt:],
                                 errors[startInt:]))
    fm = median(fluxes[startInt:])
    cm = median(counts[startInt:])
    medAndErrs.append((fm,sqrt(cm)*(fm/cm)))

    # filtering out LCs where local medians are too different from the global median
    #medCutoff = 2 # how many sigma different a local median must be to be considered too different
    #maxDiffLocMed = 0 # maximum segments that may have different local medians
    
    #difMeds = [x for x in [np.abs(med - fluxMed)/merr for med,merr in medAndErrs] if x > medCutoff]
    #if len(difMeds) > maxDiffLocMed:
    #    print("Local medians are inconsistant, ending flare detection.")
    #    return 0,[]
    
    flareLst = [] # format will be [(starTime, peakTime,endTime, peakSigma)]
    quiFLuxes = []

    for timeLst,fluxLst,errorLst in timesFluxesErrors:
        
        flarePeakInd = argmax(fluxLst)
        lcMaxInd = len(timeLst) - 1
        
        # discard interval if no flares
        if not (fluxLst[flarePeakInd] - fluxMed) > (3.5 * errorLst[flarePeakInd]):
            quiFLuxes += fluxLst
            continue

        moreFlares = True
        while(moreFlares):

            # get information about flare
            flarePeakFlux = fluxLst[flarePeakInd]
            flarePeakTime = timeLst[flarePeakInd]

            # get flare peak sigma
            flarePeakSig = (flarePeakFlux - fluxMed)/errorLst[flarePeakInd]
            
            # get flare start and end
            (flareLeftEdgeInd, flareRightEdgeInd) = findFlareEdges(fluxLst, flarePeakInd, fluxMed)

            # here I need to check if the entire interval is reading as a flare
            if (flareLeftEdgeInd == 0) and (flareRightEdgeInd == lcMaxInd):
                #print("Whole epoch was a flare, discarding epoch ")
                moreFlares = False
                continue
                

            # add flare info to flareLst
            flareLst.append((timeLst[flareLeftEdgeInd],flarePeakTime,
                             timeLst[flareRightEdgeInd],flarePeakSig))

            if flareRightEdgeInd == (len(timeLst) - 1):
                # remove flare from LC
                timeLst =   timeLst[:flareLeftEdgeInd]
                fluxLst =   fluxLst[:flareLeftEdgeInd]
                errorLst = errorLst[:flareLeftEdgeInd]
            else:
                # remove flare from LC
                timeLst =   timeLst[:flareLeftEdgeInd] +  timeLst[flareRightEdgeInd:]
                fluxLst =   fluxLst[:flareLeftEdgeInd] +  fluxLst[flareRightEdgeInd:]
                errorLst = errorLst[:flareLeftEdgeInd] + errorLst[flareRightEdgeInd:]

            # if the list is empty we're done (but there's probably a problem)
            if len(fluxLst) == 0:
                #print("Oops, the whole light curve was flares...")
                #moreFlares = False
                break
            
            # check for more flares
            flarePeakInd = argmax(fluxLst)
            if not (fluxLst[flarePeakInd] - fluxMed) > (3.5 * errorLst[flarePeakInd]):
                moreFlares = False
                quiFLuxes += fluxLst
                
    # calculate quiescent flux
    quiFlux = np.mean(quiFLuxes)

    # return info
    return quiFlux, flareLst

    
def collectFlareInfo(LCfilenames):
    if not type(LCfilenames) == list:
        LCfilenames = [LCfilenames]

    flareStats = []

    for LCdatafile in LCfilenames:
        data = genfromtxt(LCdatafile,delimiter=',',names=True)
        data = data[np.isnan(data['t_mean']) == False]
        data = data[np.isnan(data['counts']) == False]
        data = data[np.isnan(data['flux']) == False]
 
        quiFlux, flareLst = getFlareWidths(data['t_mean'].tolist(),data['flux'].tolist(),
                                           data['counts'].tolist(),1600)
        if not len(flareLst):
            os.rename(LCdatafile,galexLCPath + "nonFlaring/"+LCdatafile)
            continue
        elif quiFlux == 0:
            os.rename(LCdatafile,galexLCPath + "suspectLCs/"+LCdatafile)
            continue

        for flare in flareLst:
            flareStats.append([LCdatafile,quiFlux,flare[0],flare[2],
                               data['flux'][where(data['t_mean'] == flare[1])[0][0]],
                               median(data['flux'].tolist()),flare[3]])

    return flareStats


def peakRatioHistograms(flareStats, inline=False):
    
    peakOverMedian = [int(x[4]/x[5]) for x in flareStats]
    peakOverQuiescent = [int(x[4]/x[1]) for x in flareStats]
      
    pltTitle =  'Peak Flux Distributions'

    figFont = {'fontname':'Georgia', 'size':'16'}
    
    fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(16, 8) )
    fig.set_facecolor('w')
    fig.canvas.set_window_title(pltTitle)
    
    dbins = range(0,100,1)
    ax0.hist(peakOverMedian, dbins, facecolor='#3d8500')
    ax0.set_xlabel("Peak Flux/Median Flux", **figFont)
    ax0.set_ylabel("Number of Flares", **figFont)

    pbins = range(0,100,1) #+ [1000]
    ax1.hist(peakOverQuiescent,pbins,facecolor='#2b62d7')
    ax1.set_xlabel("Peak Flux/Quiescent Flux", **figFont)

    plt.suptitle(pltTitle, **figFont)
    
    if not inline:
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
   
    plt.show()

    return 



def flareHistograms(flareStats, inline=False):
    flareDurs = [int(x[3] - x[2]) for x in flareStats]

    sigVals = [int(x[6]) for x in flareStats]
        
    pltTitle = 'Flare Duration and Peak Flux Distributions'

    figFont = {'fontname':'Georgia', 'size':'16'}
    
    fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(16, 8))
    fig.set_facecolor('w')
    fig.canvas.set_window_title(pltTitle)
    
    dbins = range(0,800,5)
    ax0.hist(flareDurs, dbins, facecolor='#810085')
    ax0.set_xlabel("Flare Duration (s)", **figFont)
    ax0.set_ylabel("Number of Flares", **figFont)

    pbins = range(0,50,1) #+ [1000]
    ax1.hist(sigVals,pbins,facecolor='#bf0046')
    ax1.set_xlabel("Peak Flux Sigma Above Median", **figFont)
    
    plt.suptitle(pltTitle, **figFont)
    
    if not inline:
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
    
    plt.show()

    return 
    
    

def graphFlaresLC(LCfilenames, galexIds):
    if not type(LCfilenames) == list:
        LCfilenames = [LCfilenames]
    if not type(galexIds) == list:
        galexIds = [galexIds]

    if not len(galexIds) == len(LCfilenames):
        print("Length of filenams and galex ids must match.")
        return

    for LCdatafile,galexID in zip(LCfilenames,galexIds):
        
        data = genfromtxt(LCdatafile,delimiter=',',names=True)                                
 
        timeShifted = [x - data['t_mean'][0] for x in data['t_mean']]
        
        quiFlux, flareLst = getFlareWidths(timeShifted,data['flux'].tolist(),
                                           data['counts'].tolist(),1600)
        
        error = [sqrt(x)*(y/x) for x,y in zip(data['counts'],data['flux'])]

        f, ax = plt.subplots(figsize=(16, 8))
        f.set_facecolor('w')
        f.canvas.set_window_title(galexID) 

        ax.errorbar(timeShifted,data['flux'], yerr=error,fmt='-o', color='#460099')
        
        ax.axhline(median(data['flux']), color='#008d2d')
        ax.axhline(quiFlux,color='#005fdf')

        for start,peakTime,end,peakSigma in flareLst:
            ax.axvline(start,color='#ff551b')
            ax.axvline(end,color='#ff551b')
            #ax.plot(peakTime,data['flux'][where(np.array(timeShifted) == peakTime)[0][0]] ,'go')

            
            plt.suptitle("Galex ID: " + galexID,fontsize=20,fontweight='bold')
            plt.xlabel('Time (Seconds)',fontsize=20)
            plt.ylabel('Flux',fontsize=20)
        
        #mng = plt.get_current_fig_manager()
        #mng.resize(*[x//2 for x in mng.window.maxsize()])
        
        plt.show()
        
        

def readInLcstats(statfile):
    lcstats = []

    STFLE = open(statfile,'r')
    arr1 = STFLE.readlines()
    STFLE.close()
    
    sarrs = [x.split(',') for x in arr1]

    for ln in sarrs:
        lcstats.append([ln[0]]+[float(x) for x in ln[1:]])
    
    return lcstats


def saveFlareStats(flareStats, statFile):

    STFLE = open(statFile,'w')

    STFLE.write("LC Filename,Quiescent Flux, Start Time, End Time, Peak Flux, Median Flux, Peak Flux Sigma Above Median\n")
    for flare in flareStats:
        STFLE.write(','.join([str(x) for x in flare]) + '\n')

    STFLE.close()
    

def loadFlareStats(statsFile):
    return ascii.read(statsFile)


def saveDistances(distances,distFile):
    DFLE = open(distFile,'w')

    DFLE.write("Kepler ID, Distance (cm)\n")
    for ln in distances:
        DFLE.write((','.join([str(x) for x in ln])).replace('None','') + '\n')

    DFLE.close()
    
    

# Teff = effective temp in K
# Rstar = stellar radius in terms of solar radius
# return distance in m
def starDist(Teff,Rstar,Av,BCv,mv):

    Rstar = Rstar * R_sun 
    Lstar = 4 * np.pi * sigma_sb * (Rstar**2) * (Teff**4)
    MbolStar = Mbol_sun - (2.5 * log10(Lstar.value/L_sun.value))
    Mv = MbolStar - BCv
    dist = 10**(0.2*(mv - Mv + 5 - Av))

    return u.pc.to(u.cm,dist) # converting distance from parsecs to cm
    

# note limit of 10,000 for MAST query, so will need multiple files
def makeMASTInputFile():

    # not most efficient, but putting the whole file in an array
    with open(dataPath + "KGSimultaneousObservations.csv",'r') as KGOBJ:
        objStrings = KGOBJ.readlines()[1:] # first line is col headers

    objStrLst=[objStrings[x:x+10000] for x in xrange(0, len(objStrings), 10000)]

    i = 0
    for obLst in objStrLst:
        MIF = open(dataPath + 'MASTInput_KIDs_' + str(i) + '.txt','w')
        for obj in obLst:
            MIF.write(obj.split(',')[0] + '\n');
        MIF.close()
        i +=1

def getApMag(ra,dec):

    # for grepping we want the string representation of the ra to 4 decimal places
    strRa = str(round(ra,4))

    # get the lines that match the ra to 4 decimal places
    grepRes = subprocess.Popen(['grep',strRa,dataPath + "EHK2012catalog.dat"],
                               stdout=subprocess.PIPE)
    obLns = grepRes.stdout.read().split('\n')
    obLns.remove('')
    if len(obLns) == 0:
        return (None,None,None, None)

    obLsts = [x.split(',') for x in obLns]
    if len(obLsts) == 1:
        return (obLsts[0][0],obLsts[0][1],obLsts[0][6],obLsts[0][7])

    obData = obLsts[np.array([np.abs(float(x[1]) - dec) for x in obLsts]).argmin()]
    return (obData[0], obData[1], obData[6],obData[7])


# Do all of them all at once, so don't need to do again    
def matchKIDapmag():
    
    mastTable = Table.read(dataPath + 'MAST_KIC_info.csv', format='ascii')

    calcMags = []

    for kidInfo in mastTable[1:]: # first row just has data types in it
        (ra,dec, apMag, magErr) = getApMag(float(kidInfo['RA (J2000)']),
                                   float(kidInfo['Dec (J2000)']))
        calcMags.append([kidInfo['Kepler ID'],ra,dec,apMag,magErr])
        print(ra, ' ',dec, ' ', apMag, ' ', magErr)

    # writing the output
    with open(dataPath + 'EHK2012catalog_callculatedMagnitudes.csv','w') as APMAGS:
        APMAGS.write("Kepler ID,RA (J2000),DEC (J2000), Apparent Magnitude, App Mag Uncertainty\n")
        for obj in calcMags:
            APMAGS.write(','.join([str(x) for x in obj]) + '\n')
            
            
# distances in cm
def calculateDistances():

    # collecting all the info we need
    mastTable = ascii.read(dataPath + 'MAST_KIC_info.csv')
    EHKCat = ascii.read(dataPath + 'EHK2012catalog_callculatedMagnitudes.csv')
    pTable = ascii.read(dataPath + 'picaultEtAl_Table5.txt')

    missingVals = {"Teff":0,"Radius":0,"A_V":0, "m_v":0}
    parallaxStars = 0
    nonStars = []
    distances = []
    
    for objRow in mastTable:

        if objRow['Star/Gal ID'] == 1: # Object is likely galaxy, not star
            nonStars.append(objRow['Kepler ID'])
            dist = None  
        elif objRow['Parallax']: # distance is easy with paralax
            dist = u.pc.to(u.cm, 1/(objRow['Parallax'])) # in cm
            parallaxStars += 1
        elif not (objRow['Teff'] or objRow['Radius'] or objRow['A_V']):
            if not objRow['Teff']:
                missingVals['Teff'] += 1
            if not objRow['Radius']:
                missingVals['Radius'] += 1
            if not objRow['A_V']:
                missingVals['A_V'] += 1
            dist = None
        else:
            BCv = pTable['BCV'][(np.abs(pTable['Teff']-objRow['Teff'])).argmin()]
            m_v = EHKCat['Apparent Magnitude'][np.where(EHKCat['Kepler ID'] == objRow['Kepler ID'])][0]
            if not m_v:
                missingVals['m_v'] += 1
                dist = None
            else:
                dist = starDist(objRow['Teff'],objRow['Radius'],objRow['A_V'],BCv,float(m_v))

        distances.append((objRow['Kepler ID'],dist))

    print("Number of objects that are probably not stars:", len(nonStars))
    print("Number of missing values:")
    print(missingVals)
    print("Number of stars with parallax:",parallaxStars)
    return distances,nonStars


    


def makeKeplerMASTparamFile(LCfilenames):
    if not type(LCfilenames) == list:
        LCfilenames = [LCfilenames]
        
    # we are going to need to know the KIDs that go with our files
    k2g = np.genfromtxt(dataPath + "kepler2galex.csv",delimiter=',',names=True,dtype=int)

    gids = [x[:-7] for x in LCfilenames]
    kids = [str(k2g['kic_kepler_id'][np.where(k2g['glx_objid'] == int(x))][0]) for x in gids]
    
    with open("MAST_LC_input.txt",'w') as PARAMFLE:
        for k in kids:
            PARAMFLE.write(k + '\n')


def fluxVsDuration(flareStats):
    flareDurs = [x[3] - x[2] for x in flareStats]
    peakVals = [x[4]/x[1] for x in flareStats]
    fileNames = [x[0] for x in flareStats]

    figFont = {'fontname':'Georgia', 'size':'16'}
    
    fig = plt.figure(figsize=(16,10))
    ax1 = fig.add_subplot(111)
    col = ax1.loglog(flareDurs,peakVals,linestyle='None',marker='o',color='#4c1c85',alpha=.25)

    plt.xlabel('Flare Duration (seconds)',**figFont)
    plt.ylabel('Peak flux (max/quiescent flux)',**figFont)
    plt.title('Flare Duration vs Peak Flux',**figFont)
    
    af =  AnnoteFinder(flareDurs,peakVals,fileNames, ax=ax1, xtol=5,ytol=5)
    fig.canvas.mpl_connect('button_press_event', af)

    plt.show()



from scipy.fftpack import fft, fftfreq, fftshift

# NOTE: need to deal with data not always having equal time bins
# (or maybe just ignore for now.... :P
def testForPeriodicity(counts):
    N = len(counts)
    T = 1/N

    y = counts
    yf = fft(y)
    xf = fftfreq(N, T)
    xf = fftshift(xf)
    yplot = fftshift(yf)

    # we are only interested in the values > 0
    # (graph is symmetrical about 0, and there is
    # always a spike at 0)
    #initInd = np.where(xf == 0)[0][0] + 1

    plt.plot(xf, 1.0/N * np.abs(yplot))
    plt.ylim([0,200])
    plt.grid()
    plt.show()


def calculatePDM(timesCounts,minPeriod,maxPeriod,binWidth):

    varianceByPeriod = []
    
    for period in range(minPeriod,maxPeriod+1):
        countByPhase = {}
        for x in range(0,period,binWidth):
            countByPhase[x] = []
        ps = countByPhase.keys()
        for time,count in timesCounts:
            exps = time % period
            phase = ps[np.argmin([np.abs(x - exps) for x in ps])]
            countByPhase[phase] = countByPhase.get(phase,[]) + [count]
        #print(countByPhase)
        totVar = 0
        for phase in countByPhase.keys():
            if len(countByPhase[phase]) < 2:
                continue
            pmean = np.mean(countByPhase[phase])
            var = np.abs((sum([(x - pmean)**2 for x in countByPhase[phase]]))/(len(countByPhase[phase])))
            totVar += var#np.abs((sum([(x - pmean) for x in countByPhase[phase]]))/(len(countByPhase[phase]) - 1))
            #print(var)
        varianceByPeriod.append((totVar,period))
        #print(totVar,' ',period)

    return #min(varianceByPeriod),max(varianceByPeriod)
            
            


def manualFlareIdentification(lcfiles):

    if not type(lcfiles) == list:
        LCfilenames = [lcfiles]

    flareFiles = []
    notFlareFiles = []

    for lcfle in lcfiles:
        graphLightCurveFlux(lcfle)
        isFlare = raw_input("Flaring star? (y/n) ")
        if isFlare == 'y':
            flareFiles.append(lcfle)
        else:
            notFlareFiles.append(lcfle)

    return flareFiles,notFlareFiles


from gatspy.periodic import LombScargleFast

def lombScargle(LCFilename):

    data = np.genfromtxt(dataPath + 'flareFiles/' + LCFilename,delimiter=',',names=True)
    error = np.array([np.sqrt(x)*(y/x) for x,y in zip(data['counts'],data['flux'])])

    #fig, ax1 = plt.subplots()
    #ax1.errorbar(data['t_mean'], data['flux'], error, fmt='.k', ecolor='gray')
    #ax1.set(xlabel='Time (sec)', ylabel='flux',
    #       title=LCFilename)
    #ax1.invert_yaxis();
    #plt.show()

    model = LombScargleFast().fit(data['t_mean'], data['flux'], error)
    periods, power = model.periodogram_auto(nyquist_factor=100)

    fig, ax2 = plt.subplots()
    ax2.plot(periods, power)
    ax2.set(xlim=(5,1000),# ylim=(0, 0.8),
           xlabel='period (secs)',
           ylabel='Lomb-Scargle Power');

    #posPers = [(x,y) for x,y in zip(periods,power) if ((x > 5) and (x < 200))]
    #if len(posPers) == 0:
    #    period = "None"
    #else:
    #    period,power = max(posPers,key=lambda x:x[1])
    

    model.optimizer.period_range=(30, 200)
    period = model.best_period
    print("period = {0}".format(period))

    plt.show()


def saveLCFlareFigs(lcfiles):

    # turn lcfiles into list if it isn't
    if not type(LCfilenames) == list:
        LCfilenames = [LCfilenames]
        
    for LCFile in lcfiles:

        # data = etc
        data = genfromtxt(LCdatafile,delimiter=',',names=True)
        
        # error = etc
        error = [sqrt(x)*(y/x) for x,y in zip(data['counts'],data['flux'])]
        plt.errorbar(data['t_mean'],data['flux'], yerr=error,fmt='b-o')
        
        plt.axhline(median(data['flux']), color='g')
        plt.axhline(quiFlux,color='r')

        for start,peakTime,end in flareLst:
            plt.axvline(start,color='r')
            plt.axvline(end,color='r')
            plt.plot(peakTime,data['flux'][where(data['t_mean'] == peakTime)[0][0]] ,'go')

        plt.xlabel('t_mean')
        plt.ylabel('flux')
        plt.title(LCdatafile)

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()
        
        


    

def graphLightCurveFlux(LCfilenames):
    if not type(LCfilenames) == list:
        LCfilenames = [LCfilenames]

    for LCdatafile in LCfilenames:
        galexID = LCdatafile[-26:-7]

        # get full light curve
        data = genfromtxt(LCdatafile,delimiter=',',names=True)

        quiFlux, flareLst = getFlareWidths(data['t_mean'].tolist(),
                                           data['flux'].tolist(),
                                           data['counts'].tolist(),1600)

        # get galex time intervals
        obsInts = findIntervals(galexID)

        med = median(data['flux'])
        
        # make the tables for graphing
        countByTime=[]
        ii = 0
        for start,stop in obsInts:
            fi = bisect_right(data['t_mean'],stop)
            countByTime.append((data['t_mean'][ii:fi],
                                data['counts'][ii:fi],
                                data['flux'][ii:fi]))
            ii = fi

        # set up the subplots
        nplots = len(countByTime)
        if not nplots:
            print("Find intervals says there are not plots to make.")
            continue
        ncols = np.floor(sqrt(nplots)).astype(int)
        nrows = np.ceil(nplots/ncols).astype(int)

        #print(flareLst)

        f, axs = plt.subplots(nrows, ncols, figsize=(16, 8))
        f.set_facecolor('w')
        f.canvas.set_window_title(galexID)
        
        titleFont = {'fontname':'Georgia', 'size':'20'}
        axisFont = {'fontname':'Georgia', 'size':'14'}

        ax1 = None
        for i in range(nplots):
            if not ax1:
                ax1 = plt.subplot(nrows, ncols, i+1)
            else:
                plt.subplot(nrows, ncols, i+1, sharey=ax1)
            
            error = [sqrt(x)*(y/x) for x,y in zip(countByTime[i][1],countByTime[i][2])]
            plt.errorbar(countByTime[i][0], countByTime[i][2], yerr=error,fmt='-o',color='#460099')
            plt.axhline(median(countByTime[i][2]), color='#008d2d')
            #plt.axhline(med, color='r',linewidth=2)
            plt.axhline(quiFlux,color='#005fdf')

            for start,peakTime,end,peakSig in flareLst:
                if start > countByTime[i][0][-1]:
                    break
                elif end < countByTime[i][0][0]:
                    continue
                plt.axvline(start,color='#ff551b')
                plt.axvline(end,color='#ff551b')
                
            plt.xlabel('t_mean',**axisFont)
            plt.ylabel('flux (poisson error)',**axisFont)

        plt.suptitle('GALEX object id: ' + galexID,**titleFont)

        plt.show()


# NOTE: this function expects full path to LC file!!
def getFluence(LCFilepath,flareStart,flareEnd, quiescentFlux):

    data = ascii.read(LCFilepath)

    startInd = np.array([np.abs(x - flareStart) for x in data['t_mean'] if not ma.is_masked(x)]).argmin()
    endInd = np.array([np.abs(x - flareEnd) for x in data['t_mean'] if not ma.is_masked(x)]).argmin()

    intFlux = 0;
    for i in range(startInd,endInd + 1):
        if ma.is_masked(data['t_mean'][i]) or ma.is_masked(data['flux'][i]):
            continue
        if (data['flux'][i] - quiescentFlux) < 0:
            continue
        intFlux += (data['flux'][i] - quiescentFlux)*data['exptime'][i]*galexEqWth

    return intFlux



# note limit of 10,000 for MAST query, so will need multiple files
def makeMASTInputFile():

    # not most efficient, but putting the whole file in an array
    with open(dataPath + "KGSimultaneousObservations.csv",'r') as KGOBJ:
        objStrings = KGOBJ.readlines()[1:] # first line is col headers

    objStrLst=[objStrings[x:x+10000] for x in xrange(0, len(objStrings), 10000)]

    i = 0
    for obLst in objStrLst:
        MIF = open(dataPath + 'MASTInput_KIDs_' + str(i) + '.txt','w')
        for obj in obLst:
            MIF.write(obj.split(',')[0] + '\n');
        MIF.close()
        i +=1

def getApMag(ra,dec):

    # for grepping we want the string representation of the ra to 4 decimal places
    strRa = str(round(ra,4))

    # get the lines that match the ra to 4 decimal places
    grepRes = subprocess.Popen(['grep',strRa,dataPath + "EHK2012catalog.dat"],
                               stdout=subprocess.PIPE)
    obLns = grepRes.stdout.read().split('\n')
    obLns.remove('')
    if len(obLns) == 0:
        return (None,None,None, None)

    obLsts = [x.split(',') for x in obLns]
    if len(obLsts) == 1:
        return (obLsts[0][0],obLsts[0][1],obLsts[0][6],obLsts[0][7])

    obData = obLsts[np.array([np.abs(float(x[1]) - dec) for x in obLsts]).argmin()]
    return (obData[0], obData[1], obData[6],obData[7])


# Do all of them all at once, so don't need to do again    
def matchKIDapmag():
    
    mastTable = Table.read(dataPath + 'MAST_KIC_info.csv', format='ascii')

    calcMags = []

    for kidInfo in mastTable[1:]: # first row just has data types in it
        (ra,dec, apMag, magErr) = getApMag(float(kidInfo['RA (J2000)']),
                                   float(kidInfo['Dec (J2000)']))
        calcMags.append([kidInfo['Kepler ID'],ra,dec,apMag,magErr])
        print(ra, ' ',dec, ' ', apMag, ' ', magErr)

    # writing the output
    with open(dataPath + 'EHK2012catalog_callculatedMagnitudes.csv','w') as APMAGS:
        APMAGS.write("Kepler ID,RA (J2000),DEC (J2000), Apparent Magnitude, App Mag Uncertainty\n")
        for obj in calcMags:
            APMAGS.write(','.join([str(x) for x in obj]) + '\n')



# so can look at interesting ones
# returns array of the form [(LCFilename.startTime,endTime,energy),...]
def getEnergies(flareStats):

    flareEnergies = []

    # we are going to need to know the KIDs that go with our files
    k2g = ascii.read(dataPath + 'kepler2galex.csv')
        
    allDists = ascii.read(dataPath + "calculatedStarDistances.csv")
    
    for flareFle,quiFlux,start,end,peakFlux,medFlux,peakSigma in flareStats:
        gid = int(flareFle[-26:-7]) 
        kid = k2g['kic_kepler_id'][np.where(k2g['glx_objid'] == gid)][0]
        dist = allDists['Distance (cm)'][np.where(allDists['Kepler ID'] == kid)][0]

        if not dist:
            continue # can't calculate energy without distance
        fluence = getFluence(flareFle,start,end,quiFlux)
        energy =  (4 * np.pi * (dist**2) * fluence) / galexBPpercentE

        flareEnergies.append((flareFle,kid,start,end,energy))
        
    return flareEnergies


def saveEnergies(flareEnergies, filename):
    with open(filename,'w') as EFLE:
        EFLE.write("Filename, Kepler ID, Start Time, End Time, Estimated Bolometric Energy (erg)\n")
        for ln in flareEnergies:
            EFLE.write(','.join([str(x) for x in ln]) + "\n")




def graphEnergies(energies):

    energyList = np.sort(energies)[::-1]
    
    plt.loglog(energyList,range(1,len(energyList)+1), '.', color='#600b55', linestyle='None')

    plt.title("\n".join(wrap("Integrated Energy versus Number of Flares Observed in GALEX with Greater Integrated Energies (log-log)")))
    plt.xlabel("Integrated Energy")
    plt.ylabel("N(>E)")

    plt.show()


def makeKeplerMASTparamFile(LCfilenames, paramfilename = "MAST_LC_input.txt"):
    if not type(LCfilenames) == list:
        LCfilenames = [LCfilenames]
        
    # we are going to need to know the KIDs that go with our files
    k2g = np.genfromtxt(dataPath + "kepler2galex.csv",delimiter=',',names=True,dtype=int)

    gids = [x[:-7] for x in LCfilenames]
    kids = [str(k2g['kic_kepler_id'][np.where(k2g['glx_objid'] == int(x))][0]) for x in gids]
    
    with open(paramfilename,'w') as PARAMFLE:
        for k in kids:
            PARAMFLE.write(k + '\n')

    
def getD(Emin, energyList):
    energyList = sorted(energyList)    
    minInd = bisect_right(energyList,Emin)

    energyList = energyList[minInd:]
    energyList.sort(reverse = True)

    Ntot = len(energyList)
    alpha = 1 +Ntot/sum([np.log(E/Emin) for E in energyList])
    
    S_x = [(N+1)/Ntot for N in range(len(energyList))]

    P_x = [(E/Emin)**(1-alpha) for E in energyList]

    B  = max([np.abs(S_x[i] - P_x[i]) for i in range(len(S_x))])

    return B,alpha,Emin 


def fitEnergyCurve(energyList):

    # sort the energies (and copy the array so we aren't messing with the original)
    energyList = sorted(energyList)
    
    # get rid of any negative energies (where do they come from!? argh)
    minInd = bisect_right(energyList,0)
    energyList = energyList[minInd:]

    BalphaEmin = []

    for Emin in energyList[:len(energyList)//2]:
        BalphaEmin.append(getD(Emin,energyList))

    # our result should be when B is minimized
    B,alpha,Emin = min(BalphaEmin,key=lambda tup: tup[0]) 

    print("B",B)
    print("Emin",Emin)
    print("alpha",alpha)

    return B,alpha,Emin
    

def graphFitAndEnergies(Emin,alpha,energies):

    energyList= np.sort(energies)[::-1]
    
    Ntot = len([E for E in energyList if E >= Emin])

    figFont = {'fontname':'Georgia', 'size':'16'}
    
    fig = plt.figure(figsize=(16,10))
    ax1 = fig.add_subplot(111)

    at = AnchoredText(r'$\alpha = %.3f$' '\n' r'$E_{min} = %.3E$' % (alpha, Emin),
                      prop=dict(size=18), frameon=False,
                      loc=1)
    ax1.add_artist(at)
    
    ax1.loglog(energyList,range(1,len(energyList)+1), '.', color='#600b55', linestyle='None')
    ax1.loglog(energyList,[Ntot*((E/Emin)**(1-alpha)) for E in energyList], color='#2f6108',linewidth=2)

    plt.title("\n".join(wrap("Integrated Energy versus Number of Flares Observed in GALEX with Greater Integrated Energies (log-log)")),**figFont)
    plt.xlabel("Integrated Energy",**figFont)
    plt.ylabel("N(>E)",**figFont)
    
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)   
        
    plt.show()


from scipy.fftpack import fft, fftfreq, fftshift

# NOTE: need to deal with data not always having equal time bins
# (or maybe just ignore for now.... :P
def testForPeriodicity(counts):
    N = len(counts)
    T = 1/N

    y = counts
    yf = fft(y)
    xf = fftfreq(N, T)
    xf = fftshift(xf)
    yplot = fftshift(yf)

    # we are only interested in the values > 0
    # (graph is symmetrical about 0, and there is
    # always a spike at 0)
    #initInd = np.where(xf == 0)[0][0] + 1

    plt.plot(xf, 1.0/N * np.abs(yplot))
    plt.ylim([0,200])
    plt.grid()
    plt.show()


def calculatePDM(timesCounts,minPeriod,maxPeriod,binWidth):

    varianceByPeriod = []
    
    for period in range(minPeriod,maxPeriod+1):
        countByPhase = {}
        for x in range(0,period,binWidth):
            countByPhase[x] = []
        ps = countByPhase.keys()
        for time,count in timesCounts:
            exps = time % period
            phase = ps[np.argmin([np.abs(x - exps) for x in ps])]
            countByPhase[phase] = countByPhase.get(phase,[]) + [count]
        #print(countByPhase)
        totVar = 0
        for phase in countByPhase.keys():
            if len(countByPhase[phase]) < 2:
                continue
            pmean = np.mean(countByPhase[phase])
            var = np.abs((sum([(x - pmean)**2 for x in countByPhase[phase]]))/(len(countByPhase[phase])))
            totVar += var#np.abs((sum([(x - pmean) for x in countByPhase[phase]]))/(len(countByPhase[phase]) - 1))
            #print(var)
        varianceByPeriod.append((totVar,period))
        #print(totVar,' ',period)

    return min(varianceByPeriod),max(varianceByPeriod)
            
            


#def manualFlareIdentification(lcfiles):

#    if not type(lcfiles) == list:
#        LCfilenames = [lcfiles]

#    flareFiles = []
#    notFlareFiles = []

#    for lcfle in lcfiles:
#        graphLightCurveFlux(lcfle)
#        isFlare = raw_input("Flaring star? (y/n) ")
#        if isFlare == 'y':
#            flareFiles.append(lcfle)
#        else:
#            notFlareFiles.append(lcfle)

#    return flareFiles,notFlareFiles


from gatspy.periodic import LombScargleFast

def lombScargle(LCFilename):

    data = np.genfromtxt(dataPath + 'flareFiles/' + LCFilename,delimiter=',',names=True)
    error = np.array([np.sqrt(x)*(y/x) for x,y in zip(data['counts'],data['flux'])])

    #fig, ax1 = plt.subplots()
    #ax1.errorbar(data['t_mean'], data['flux'], error, fmt='.k', ecolor='gray')
    #ax1.set(xlabel='Time (sec)', ylabel='flux',
    #       title=LCFilename)
    #ax1.invert_yaxis();
    #plt.show()

    model = LombScargleFast().fit(data['t_mean'], data['flux'], error)
    periods, power = model.periodogram_auto(nyquist_factor=100)

    fig, ax2 = plt.subplots()
    ax2.plot(periods, power)
    ax2.set(xlim=(5,1000),# ylim=(0, 0.8),
           xlabel='period (secs)',
           ylabel='Lomb-Scargle Power');

    #posPers = [(x,y) for x,y in zip(periods,power) if ((x > 5) and (x < 200))]
    #if len(posPers) == 0:
    #    period = "None"
    #else:
    #    period,power = max(posPers,key=lambda x:x[1])
    

    model.optimizer.period_range=(30, 200)
    period = model.best_period
    print("period = {0}".format(period))

    plt.show()


def saveLCFlareFigs(LCfilenames):

    # turn lcfiles into list if it isn't
    if not type(LCfilenames) == list:
        LCfilenames = [LCfilenames]
        
    for LCFile in LCfilenames:
        #print(LCFile)
        
        # data = etc
        data = genfromtxt(LCFile,delimiter=',',names=True)
        med = median(data['flux'])
        
        # error = etc
        error = [sqrt(x)*(y/x) for x,y in zip(data['counts'],data['flux'])]
        
        # get the flares
        quiFlux, flareLst = getFlareWidths(data['t_mean'].tolist(),
                                           data['flux'].tolist(),
                                           data['counts'].tolist(),1600)
        
        # if there aren't any flares onto the next
        if len(flareLst) == 0:
            continue
        
        #print(len(flareLst))
        
        # get the intervals
        obsInts = findIntervals(LCFile[:-7])

        # sepparate the intervals
        fluxByTime=[]
        ii = 0
        for start,stop in obsInts:
            fi = bisect_right(data['t_mean'],stop)
            fluxByTime.append((data['t_mean'][ii:fi],
                               data['flux'][ii:fi],
                               error[ii:fi]))
            ii = fi

        #print(len(obsInts))
        # for each flare:
        for startTime,peakTime,endTime in flareLst:
            
            # determine which interval it is in
            for times,fluxes,errors in fluxByTime:
                if endTime <= times[-1]:
                    break
                    
            # set up graph of just that interval, with single flare markings
            plt.clf()
            plt.errorbar(times, fluxes, yerr=errors,fmt='b-o')
            plt.axhline(median(fluxes), color='g',linewidth=2)
            plt.axhline(med, color='r',linewidth=2)
            plt.axhline(quiFlux,color='b',linewidth=2)
            plt.axvline(startTime,color='r')
            plt.axvline(endTime,color='r')
            #plt.show()
        
            # instead of displaying save as:
            plt.savefig(dataPath + "flareFiles/flarePngs/" + LCFile[:-7] + "_" + str(startTime) + "_" + str(endTime) + ".png")
            
        

def getGoodFlareStats(pngFilenames,masterStatList):

    pngStats = [x[:-4].split('_') for x in pngFilenames]

    goodStats = []

    for goodPng in pngStats:
        startTime = goodPng[1]
        for stat in masterStatList:
            if str(stat[2]) == startTime:
                goodStats.append(stat)
                break

    return goodStats

from astropy.time import Time
def graphFlare(flareStat, galexID = None):
    if not galexID:
        galexID = flareStat[0][-26:-7]
        
    data = genfromtxt(flareStat[0],delimiter=',',names=True)
    med = median(data['flux'])
        
    error = [sqrt(x)*(y/x) for x,y in zip(data['counts'],data['flux'])]
        
    # get the intervals
    obsInts = findIntervals(galexID)

    # separate the intervals
    fluxByTime=[]
    ii = 0
    for start,stop in obsInts:
        fi = bisect_right(data['t_mean'],stop)
        fluxByTime.append((data['t_mean'][ii:fi],
                               data['flux'][ii:fi],
                               error[ii:fi]))
        ii = fi

    startTime = flareStat[2]
    endTime = flareStat[3]
    quiFlux = flareStat[1]

    
    # determine which interval it is in
    for times,fluxes,errors in fluxByTime:
        if len(times) == 0:
            continue
        if endTime <= times[-1]:
            break
                    
    # set up graph of just that interval, with single flare markings
    plt.clf()
    plt.errorbar(times, fluxes, yerr=errors,fmt='b-o')
    plt.axhline(median(fluxes), color='g',linewidth=1)
    plt.axhline(med, color='r',linewidth=1)
    plt.axhline(quiFlux,color='b',linewidth=2)
    plt.axvline(startTime,color='r')
    plt.axvline(endTime,color='r')
    plt.show()
        
          
def getAllInfoForFlare(LCFile,flareStart,flareEnd):
    pass

def graphFitsLCs(fitsLCFiles):

    if not type(fitsLCFiles) == list:
        fitsLCFiles = [fitsLCFiles]
        
    for fitsFile in fitsLCFiles:
        kdata = fits.open(fitsFile)

        lcTab = kdata['LIGHTCURVE'].data

        plt.clf()
        plt.errorbar(lcTab['TIME'], lcTab['SAP_FLUX'], yerr=lcTab['SAP_FLUX_ERR'],fmt='b-o')
    
        plt.title(fitsFile)
        plt.xlabel("Time (Kepler Barycentric Julian Day)")
        plt.ylabel("Flux (electrons/sec)")

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        
        plt.show()


from astropy.time import Time

# TODO: allow for multiple fits files
def graphGalexAndKepler_draft1(galexLCFile,KeplerFitsFile):

    # get Galex data
    galexData = genfromtxt(galexLCFile,delimiter=',',names=True)
    gErrors = [sqrt(x)*(y/x) for x,y in zip(galexData['counts'],galexData['flux'])]

    # get Kepler data
    fitsFile = fits.open(KeplerFitsFile)
    keplerData = fitsFile['LIGHTCURVE'].data

    # get Galex intervals
    obsInts = findIntervals(galexLCFile[:-7])

    # translate Galex intervals into astropy time
    # "GALEX Time" = "UNIX Time" - 315964800
    obsIntStdTime = []
    for start,stop in obsInts:
        obsIntStdTime.append([Time(start+315964800,format='unix'),
                              Time(stop+315964800,format='unix')])

 
    # put data all together
    galexTimeFluxError=[]
    keplerTimeFluxError=[]
    ii = 0
    for start,stop in obsInts:
        fi = bisect_right(galexData['t_mean'],stop)
        galexTimeFluxError.append((galexData['t_mean'][ii:fi],
                            galexData['flux'][ii:fi],
                            gErrors[ii:fi]))
        ii = fi

    ii = obsIntStdTime[0][0].jd
    for start,stop in obsIntStdTime:
        fi = bisect_right(keplerData['TIME'],stop.jd)
        keplerTimeFluxError.append((keplerData['TIME'][ii:fi],
                            keplerData['SAP_FLUX'][ii:fi],
                            keplerData['SAP_FLUX_ERR'][ii:fi]))
        ii = fi
        
    # graph
    for i in range(len(galexTimeFluxError)):
        
        #ax1 =
        plt.subplot(2, 1, 0)
        plt.errorbar([Time(x + 315964800,format='unix').jd for x in galexTimeFluxError[i][0]],
                     galexTimeFluxError[i][1],
                     yerr=galexTimeFluxError[i][2],fmt='b-o')
            
        plt.subplot(2, 1, 1)#, sharex=ax1)
        print(len(keplerTimeFluxError[i][0]))
        plt.errorbar(keplerTimeFluxError[i][0],
                     keplerTimeFluxError[i][1],
                     yerr=keplerTimeFluxError[i][2],fmt='g-o')            
                
        plt.xlabel('Time (Julian Day)')
        plt.ylabel('Flux')

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())

        plt.show()


keplerQuarters = [Time(54953.0283,format='mjd'),
                  Time(54964.0011,format='mjd'),
                  Time(55002.0076,format='mjd'),
                  Time(55092.7123,format='mjd'),
                  Time(55184.8679,format='mjd'),
                  Time(55275.9813,format='mjd'),
                  Time(55371.9375,format='mjd'),
                  Time(55462.6626,format='mjd'),
                  Time(55567.8548,format='mjd'),
                  Time(55641.0071,format='mjd'),
                  Time(55739.3336,format='mjd'),
                  Time(55833.6959,format='mjd'),
                  Time(55931.8998,format='mjd'),
                  Time(56015.2280,format='mjd'),
                  Time(56106.6275,format='mjd'),
                  Time(56205.9756,format='mjd'),
                  Time(56304.5882,format='mjd'),
                  Time(56391.7170,format='mjd'),
                  Time(54964.0011,format='mjd')]

                      
def makeKeplerParamFiles(flareStats):
    kidsByQuarter = {}

    # we are going to need to know the KIDs that go with our files
    k2g = np.genfromtxt(dataPath + "kepler2galex.csv",delimiter=',',names=True,dtype=int)

    # collecting the KIDs we need to download by quarter
    for flare in flareStats:
        gid = int(flare[0][:-7])
        kid = str(k2g['kic_kepler_id'][np.where(k2g['glx_objid'] == gid)][0])


        flarestart = Time(flare[2] + 315964800,format='unix')
        flareend = Time(flare[3] + 315964800,format='unix')
        quarter = np.array([np.abs(x - flarestart) for x in keplerQuarters]).argmin()
        if flarestart < keplerQuarters[quarter]:
            quarter -= 1
        if (quarter == 18) and (flarestart > keplerQuarters[quarter]):
            print("No kepler data for flare time frame.")
            continue
        if (quarter == -1) and (flareend <= keplerQuarters[0]):
            print("No kepler data for flare time frame.")
            continue
        else:
            quarter += 1

        while(flareend > keplerQuarters[quarter]):
            if quarter not in kidsByQuarter:
                kidsByQuarter[quarter] = [kid]
            else:
                if kid not in kidsByQuarter[quarter]:
                    kidsByQuarter[quarter].append(kid)
            quarter += 1

            
    # making mast param files (max 10,000 lines, not dealing with this, unliekly to come up
    for quarter in kidsByQuarter:
        filename = "mast_kids_quarter_" + str(quarter) + ".param"
        MIF = open(dataPath + 'flareFiles/' + filename,'w')
        for k in kidsByQuarter[quarter]:
            MIF.write(k + '\n');
        MIF.close()
        
    

prop = font_manager.FontProperties(fname=fontpath)
matplotlib.rcParams['font.family'] = prop.get_name()


class kgGraphData:
    """Object to hold all the information necessary to graph GALEX and Kepler LCs together.
    To initialize neet to provide a start time and stop time with type astropy Time"""
    
    def __init__(self,st, et):
        self.startTime = st
        self.endTime = et

        self.galexID = None
        self.keplerID = None
        
        self.keplerData = []
        self.galexData = []

        self.starData = ""
        self.quiescentFlux = None
        self.flareLocs = []


    def addKeplerData(self,keplerID):
        self.keplerID = keplerID
        
        # get all kepler fits files for this object
        kepfiles = [x for x in os.listdir(keplerLCPath) if keplerID in x]
        if len(kepfiles) == 0:
            print("No kepler files found")
            return
    
        # get Kepler data
        for kfle in kepfiles:
            fitsFile = fits.open(keplerLCPath + kfle)
            keplerData = fitsFile['LIGHTCURVE'].data
            #kStdTimes = [Time(x + 2454833, format='jd') if not np.isnan(x) else x for x in keplerData['TIME']]
            kStdTimes = [Time(x + 2454833, format='jd') for x in keplerData['TIME'] if not np.isnan(x)]
            fluxes = [keplerData['SAP_FLUX'][i] for i in xrange(len(keplerData['TIME'])) if not np.isnan(keplerData['TIME'][i])]
            errs = [keplerData['SAP_FLUX_ERR'][i] for i in xrange(len(keplerData['TIME'])) if not np.isnan(keplerData['TIME'][i])]
            
            # checking for data in the interval we care about
            if (self.startTime > kStdTimes[-1]) or (self.endTime < kStdTimes[0]):
                continue

            # adding the starData if we need to
            if not self.starData:
                starInfo = fitsFile[0].header
                self.starData += "Right ascension (deg): " + str(starInfo['RA_OBJ'])
                self.starData += "\nDeclination (deg): " + str(starInfo['DEC_OBJ'])
                self.starData += "\nEffective temperature (K): " + str(starInfo['TEFF'])
                self.starData += "\nStellar radius (solar radii): " + str(starInfo['RADIUS'])
                self.starData += "\nKepler magnitude (mag): " + str(starInfo['KEPMAG'])

            # collecting actual Kepler data
            ii = bisect_left(kStdTimes,self.startTime)
            fi = bisect_right(kStdTimes,self.endTime)
            self.keplerData.append((kStdTimes[ii:fi],
                                    fluxes[ii:fi],
                                    errs[ii:fi]))

            
    def addGalexData(self,galexID):
        self.galexID = galexID
        
        # get Galex data
        galexData = genfromtxt(galexLCPath + galexID + '_LC.csv',delimiter=',',names=True)
        gErrors = [sqrt(x)*(y/x) for x,y in zip(galexData['counts'],galexData['flux'])]
        gStdTimes = [Time(x+315964800,format='unix') for x in galexData['t_mean']]
       
        # get Galex intervals
        obsInts = findIntervals(galexID)
        for intval in obsInts:
            start = Time(intval[0]+315964800,format='unix')
            stop = Time(intval[1]+315964800,format='unix')
            
            if (start < self.startTime) or (stop > self.endTime): # interval is not in graphData time range
                continue
                
            ii = bisect_left(gStdTimes,start)
            fi = bisect_right(gStdTimes,stop)
            self.galexData.append((gStdTimes[ii:fi],
                                   galexData['flux'][ii:fi],
                                   gErrors[ii:fi]))
        

    def addFlareLocs(self, gflares):
        for lcfilename,quiFlux,startGtime,stopGtime,peakFlux in gflares:
            start = Time(float(startGtime)+315964800,format='unix')
            stop = Time(float(stopGtime)+315964800,format='unix')
            if self.galexID not in lcfilename:
                print("Oops, flare is not on this object.")
                continue
            if (start > self.endTime) or (stop < self.startTime): # flare is not in time interval
                continue
                
            if not self.quiescentFlux:
                self.quiescentFlux = quiFlux

            # adding flare info
            self.flareLocs.append((start,stop))




def isSubInterval(subInt,superInts):
    if type(superInts) == tuple:
        superInts = [superInts]

    subStart = float(subInt[0])
    subEnd= float(subInt[1])
    #print("Subinterval start and end")
    #print(subStart,subEnd)
    #print()

    for iStart,iEnd in superInts:
        #print("Superinterval start and end")
        #print(iStart,iEnd)
        #print()

        if (subStart >= iStart) and (subEnd <= iEnd):
            #print("Found a super-interval!")
            return (iStart,iEnd)

    return False


def getBigInts(galexIntervals, galexFlares, extKepTm):

    bigTimeInts = []

    for flare in galexFlares:
        supInt = isSubInterval((flare[2],flare[3]),galexIntervals)

        if(supInt):
            newIntSt = supInt[0] - extKepTm/2
            newIntEd = supInt[1] + extKepTm/2

            # check if should add to existing big int
            bigIntStart = isSubInterval((newIntSt,newIntSt),bigTimeInts)
            bigIntEnd = isSubInterval((newIntEd,newIntEd),bigTimeInts)
            if (bigIntStart and bigIntEnd) and (bigIntStart != bigIntEnd):                
                ind1 = bigTimeInts.index(bigIntStart)
                ind2 = bigTimeInts.index(bigIntEnd)
                ns = bigTimeInts[ind1][0]
                ne = bigTimeInts[ind1][1]
                bigTimeInts.pop(ind1)
                bigTimeInts.pop(ind2)
                bigTimeInts.append((ns,ne))
                bigTimeInts.sort()
            elif bigIntStart:
                ind = bigTimeInts.index(bigIntStart)
                st = bigTimeInts[ind][0]
                bigTimeInts[ind] = (st,newIntEd)
            elif bigIntEnd:
                ind = bigTimeInts.index(bigIntEnd)
                end = bigTimeInts[ind][1]
                bigTimeInts[ind] = (newIntSt,end)
            else: # otherwise make new big interval
                bigTimeInts.append((newIntSt,newIntEd))
                bigTimeInts.sort()

    # translate Galex time into astropy time
    # "GALEX Time" = "UNIX Time" - 315964800
    for i in range(len(bigTimeInts)):
        start = Time(bigTimeInts[i][0]+315964800,format='unix')
        stop = Time(bigTimeInts[i][1]+315964800,format='unix')
        bigTimeInts[i] = (start,stop)
 
    return bigTimeInts




        
def graphGalexAndKepler_bigInts(galexID, extraKepHours = 8, markFlares = True):
    print(galexID)
    
    # get Kepler ID
    grepRes = subprocess.Popen(['grep',galexID,dataPath + "kepler2galex.csv"],stdout=subprocess.PIPE)
    keplerID = grepRes.stdout.read().split(',')[0]
    print(keplerID)

    # get Galex intervals
    obsInts = findIntervals(galexID)
    
    # get flares 
    grepRes = subprocess.Popen(['grep',galexID,dataPath + "flareFiles/vettedStats_20160304.csv"],stdout=subprocess.PIPE)
    gflares = [x.split(',') for x in grepRes.stdout.read().strip('\n').split('\n')]
    print("Number of flares:",len(gflares))
    
    # make big intervals
    gExtraTime = extraKepHours * 60 * 60 # galex time is in seconds
    bigInts = getBigInts(obsInts, gflares, gExtraTime) # bigInts are in astropy Time
    if len(bigInts) == 0:
        print("No flaring intervals... Goodbye.")
        return

    dataForGraphs = []
    for start,stop in bigInts:
        newGraphDat = kgGraphData(start,stop)
        newGraphDat.addGalexData(galexID)
        newGraphDat.addKeplerData(keplerID)
        newGraphDat.addFlareLocs(gflares)
        dataForGraphs.append(newGraphDat)

   
    # TODO:  (incorporate into object?)(not yet)
    # graph (NOTE: MAKE BETTER)
    for graphObj in dataForGraphs:
        f, axs = plt.subplots(2, sharex = True)
        f.set_facecolor('w')
        f.canvas.set_window_title(galexID + '_' + keplerID) 
        
        at = AnchoredText(graphObj.starData,
                          prop=dict(size=18), frameon=False,
                          loc=1)
        axs[0].add_artist(at)
        
        ## Galex ##
        
        # graphing the data
        label = True
        for times,data,err in graphObj.galexData:
            if label:
                axs[0].errorbar([x.jd for x in times],data,yerr=err,fmt='b-o',label='GALEX')
                label = False
            else:
                axs[0].errorbar([x.jd for x in times],data,yerr=err,fmt='b-o')

        # marking the flares
        if markFlares:
            for fmin,fmax in graphObj.flareLocs:
                axs[0].axvline(fmin.jd,color='r')
                axs[0].axvline(fmax.jd,color='r')
            
        axs[0].axhline(graphObj.quiescentFlux,color='g')
        axs[0].legend(loc='upper left')
        
        #miny = min(data)
        #maxy = max(data)
        #plt.yticks(np.arange(miny,maxy,(maxy-miny)/6))

        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        
        ## Kepler ##
        # graphing the data
        lclabel = True
        sclabel = True
        for times,data,err in graphObj.keplerData:
            if times[1].jd - times[0].jd < 0.02: # short cadence
                if sclabel:
                    axs[1].errorbar([x.jd for x in times],data,yerr=err,fmt='b-o',label='Kepler Short Cadence')
                    sclabel = False
                else:
                    axs[1].errorbar([x.jd for x in times],data,yerr=err,fmt='b-o')
            else: # long cadence
                if lclabel:
                    axs[1].errorbar([x.jd for x in times],data,yerr=err,fmt='g-o',label='Kepler Long Cadence')
                    sclabel = False
                else:
                    axs[1].errorbar([x.jd for x in times],data,yerr=err,fmt='g-o')
            #for y in (np.array([x.jd for x in times[1:]]) + np.array([x.jd for x in times[:-1]]))/2:
            #    axs[0].axvline(y,color='g',linestyle='--')
            #    axs[1].axvline(y,color='g',linestyle='--')

        
                    
        if markFlares:       
            # marking the flares
            for fmin,fmax in graphObj.flareLocs:
                axs[1].axvline(fmin.jd,color='r')
                axs[1].axvline(fmax.jd,color='r')
            
        axs[1].legend(loc='upper left')
        
        #miny = min(data)
        #maxy = max(data)
        #plt.yticks(np.arange(miny,maxy,(maxy-miny)/6))
                
        plt.xlabel('Time (Julian Day)',fontsize=20)
        plt.suptitle('Kepler ' + keplerID,fontsize=20,fontweight='bold')
        f.text(0.1, 0.5, 'Flux', ha='center', va='center', rotation='vertical',fontsize=20)

        mng = plt.get_current_fig_manager()
        mng.resize(*[x//2 for x in mng.window.maxsize()])
        
        plt.show()



def graphGalexAndKepler_inline(galexID, extraKepHours = 8, markFlares = True, markKeplerBins = False):
    
    # get Kepler ID
    grepRes = subprocess.Popen(['grep',galexID,dataPath + "kepler2galex.csv"],stdout=subprocess.PIPE)
    keplerID = grepRes.stdout.read().split(',')[0]

    # get Galex intervals
    obsInts = findIntervals(galexID)
    
    # get flares 
    grepRes = subprocess.Popen(['grep',galexID,dataPath + "flareFiles/vettedStats_20160304.csv"],stdout=subprocess.PIPE)
    gflares = [x.split(',') for x in grepRes.stdout.read().strip('\n').split('\n')]
    
    # make big intervals
    gExtraTime = extraKepHours * 60 * 60 # galex time is in seconds
    bigInts = getBigInts(obsInts, gflares, gExtraTime) # bigInts are in astropy Time
    numInts = len(bigInts)
    if numInts == 0:
        print("No flaring intervals... Goodbye.")
        return

    if numInts == 1:
        interval = 0
    else: # NOTE: add a check to make sure the imput is valid
        interval = int(raw_input("There are %d intervals, which would you like to see? " % numInts)) - 1

    start,stop = bigInts[interval]
    graphObj = kgGraphData(start,stop)
    graphObj.addGalexData(galexID)
    graphObj.addKeplerData(keplerID)
    graphObj.addFlareLocs(gflares)
   
    # TODO:  (incorporate into object?)(not yet)
    # graph (NOTE: MAKE BETTER)
    #for graphObj in dataForGraphs:
    f, axs = plt.subplots(2, sharex = True, figsize=(16,10))
    f.set_facecolor('w')
    f.canvas.set_window_title(galexID + '_' + keplerID) 
        
    at = AnchoredText(graphObj.starData,
                      prop=dict(size=18), frameon=False,
                      loc=1)
    axs[0].add_artist(at)
        
    ## Galex ##
    
    # graphing the data
    label = True
    for times,data,err in graphObj.galexData:
        if label:
            axs[0].errorbar([x.jd for x in times],data,yerr=err,fmt='-o', color='#460099',label='GALEX')
            label = False
        else:
            axs[0].errorbar([x.jd for x in times],data,yerr=err,fmt='-o', color='#460099')

    # marking the flares
    if markFlares:
        for fmin,fmax in graphObj.flareLocs:
            axs[0].axvline(fmin.jd,color='#ff551b')
            axs[0].axvline(fmax.jd,color='#ff551b')
            
    axs[0].axhline(graphObj.quiescentFlux,color='#005fdf')
    axs[0].legend(loc='upper left')
        

    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        
    ## Kepler ##
    # graphing the data
    lclabel = True
    sclabel = True
    for times,data,err in graphObj.keplerData:
        if times[1].jd - times[0].jd < 0.02: # short cadence
            if sclabel:
                axs[1].errorbar([x.jd for x in times],data,yerr=err,fmt='-o',color='#008716',label='Kepler Short Cadence')
                sclabel = False
            else:
                axs[1].errorbar([x.jd for x in times],data,yerr=err,fmt='-o',color='#008716')
        else: # long cadence
            if lclabel:
                axs[1].errorbar([x.jd for x in times],data,yerr=err,fmt='-o',color='#008716',label='Kepler Long Cadence')
                sclabel = False
            else:
                axs[1].errorbar([x.jd for x in times],data,yerr=err,fmt='g-o')
        if markKeplerBins:
            for y in (np.array([x.jd for x in times[1:]]) + np.array([x.jd for x in times[:-1]]))/2:
                axs[0].axvline(y,color='#008716',linestyle='--')
                axs[1].axvline(y,color='#008716',linestyle='--')

        
                    
    if markFlares:       
        # marking the flares
        for fmin,fmax in graphObj.flareLocs:
            axs[1].axvline(fmin.jd,color='#ff551b')
            axs[1].axvline(fmax.jd,color='#ff551b')
            
    axs[1].legend(loc='upper left')
                
    plt.xlabel('Time (Julian Day)',fontsize=20)
    plt.suptitle('Kepler ' + keplerID,fontsize=20,fontweight='bold')
    f.text(0.1, 0.5, 'Flux', ha='center', va='center', rotation='vertical',fontsize=20)
        
    plt.show()



        
def graphGalexAndKepler(galexID, extraKepHours = 8, allIntervals = 'no'):
    print(galexID)
    
    # get Galex data
    galexData = genfromtxt(galexLCPath + galexID + '_LC.csv',delimiter=',',names=True)
    gErrors = [sqrt(x)*(y/x) for x,y in zip(galexData['counts'],galexData['flux'])]

    # get Galex intervals
    obsInts = findIntervals(galexID)

    # get flares 
    grepRes = subprocess.Popen(['grep',galexID,dataPath + "flareFiles/vettedStats_20160304.csv"],stdout=subprocess.PIPE)
    gflares = [x.split(',') for x in grepRes.stdout.read().strip('\n').split('\n')]
    print("Number of flares:",len(gflares))
    
    quiFlux = float(gflares[0][1])

    if allIntervals == 'no':
        flareInts = []
        for flare in gflares:
            for itv in obsInts:
                if (float(flare[2]) >= itv[0]) and (float(flare[3]) <= itv[1]) and (itv not in flareInts):
                    flareInts.append(itv)
                    break
            if len(flareInts) == len(obsInts):
                break
    else:
        flareInts = obsInts

            
    # get galex data
    galexTimeFluxError=[]
    intNum = 0
    for start,stop in flareInts:
        ii = bisect_left(galexData['t_mean'],start)
        fi = bisect_right(galexData['t_mean'],stop)
        galexTimeFluxError.append((intNum,
                                   galexData['t_mean'][ii:fi],
                                   galexData['flux'][ii:fi],
                                   gErrors[ii:fi]))
        intNum += 1

        
    # translate Galex intervals into astropy time
    # "GALEX Time" = "UNIX Time" - 315964800
    flareIntStdTime = []
    for start,stop in flareInts:
        flareIntStdTime.append([Time(start+315964800,format='unix'),
                                Time(stop+315964800,format='unix')])

    # get Kepler ID
    grepRes = subprocess.Popen(['grep',galexID,dataPath + "kepler2galex.csv"],stdout=subprocess.PIPE)
    keplerID = grepRes.stdout.read().split(',')[0]
    print(keplerID)

    # get all kepler fits files for this object
    kepfiles = [x for x in os.listdir(keplerLCPath) if keplerID in x]
    if len(kepfiles) == 0:
        print("No kepler files found")
        return
    
    # get Kepler data
    keplerTimeFluxError = []
    starData = ""
    for kfle in kepfiles:
        fitsFile = fits.open(keplerLCPath + kfle)
        keplerData = fitsFile['LIGHTCURVE'].data
        
        if not starData:
            starInfo = fitsFile[0].header
            starData += "Right ascension (deg): " + str(starInfo['RA_OBJ'])
            starData += "\nDeclination (deg): " + str(starInfo['DEC_OBJ'])
            starData += "\nEffective temperature (K): " + str(starInfo['TEFF'])
            starData += "\nStellar radius (solar radii): " + str(starInfo['RADIUS'])
            starData += "\nKepler magnitude (mag): " + str(starInfo['KEPMAG'])
            
        flareNum = -1
        for start,stop in flareIntStdTime:
            flareNum += 1
            if (keplerData['TIME'][-1]+2454833 <= start.jd) or (keplerData['TIME'][0]+2454833 >= stop.jd):
                continue

            # 0.2 jd ~ 4 hours
            ii = bisect_left(keplerData['TIME'],start.jd - 2454833 - .05*extraKepHours/2) 
            if ii < 0:
                ii = 0
            fi = bisect_right(keplerData['TIME'],stop.jd - 2454833 + .05*extraKepHours/2)
            if fi > len(keplerData['TIME']):
                fi = len(keplerData['TIME'])
            keplerTimeFluxError.append((flareNum,
                                        keplerData['TIME'][ii:fi],
                                        keplerData['SAP_FLUX'][ii:fi],
                                        keplerData['SAP_FLUX_ERR'][ii:fi]))

    # graph (NOTE: MAKE BETTER)
    for i in range(len(flareInts)):
        
        plots = [x[1:] for x in galexTimeFluxError if x[0] == i]
        plots += [x[1:] for x in keplerTimeFluxError if x[0] == i]

        if len(plots) == 1:
            print("No Kepler data for interval, skipping")
            continue
        else:
            startTime = plots[1][0][0] + 2454833

        #plt.locator_params(axis='y',nbins=6)

        f, axs = plt.subplots(len(plots), sharex = True)
        f.set_facecolor('w')
        #f.set_figwidth(12)
        #f.set_figheight(6)
        #f.set_label("TEST")
        f.canvas.set_window_title(galexID + '_' + keplerID) 
        f.text(.75, .81, starData, fontsize=18)
        #f.text(.13, .84, "GALEX",color='b',fontsize=24,fontweight='bold')
        #f.text(.13, .82, "Kepler",color='g',fontsize=24,fontweight='bold')
        

        axs[0].errorbar([Time(x + 315964800,format='unix').jd  - startTime for x in plots[0][0]],
                        plots[0][1],yerr=plots[0][2],fmt='b-o',label='GALEX')
        axs[0].axhline(quiFlux,color='g')
        axs[0].legend(loc='upper left')
        miny = min(plots[0][1])
        maxy = max(plots[0][1])
        plt.yticks(np.arange(miny,maxy,(maxy-miny)/6))

        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        #plt.setp([a.get_xticks() for a in f.axes[:-1]], visible=False)
        
        for flare in gflares:
            if (float(flare[2]) >=  plots[0][0][0]) and (float(flare[3]) <= plots[0][0][-1]):
                axs[0].axvline(Time(float(flare[2]) + 315964800,format='unix').jd - startTime,color='r')
                axs[0].axvline(Time(float(flare[3]) + 315964800 ,format='unix').jd - startTime,color='r')
        
        for j in range(1,len(plots)):
            #plt.subplot(len(plots), 1, j+1, sharex=ax1)
            if (plots[j][0][1] - plots[j][0][0]) < 0.02:
                pltLabel = 'Kepler short cadence'
            else:
                pltLabel = 'Kepler long cadence'
            axs[j].errorbar([x + 2454833 - startTime for x in plots[j][0]],
                            plots[j][1],yerr=plots[j][2],fmt='g-o',label=pltLabel)
            axs[j].legend(loc='upper left')
            for flare in gflares:
                if (float(flare[2]) >=  plots[0][0][0]) and (float(flare[3]) <= plots[0][0][-1]):
                    axs[j].axvline(Time(float(flare[2]) + 315964800,format='unix').jd - startTime,color='r')
                    axs[j].axvline(Time(float(flare[3]) + 315964800,format='unix').jd - startTime,color='r')
            miny = min(plots[j][1])
            maxy = max(plots[j][1])
            plt.yticks(np.arange(miny,maxy,(maxy-miny)/6))

                
        plt.xlabel('Time (Julian Day)',fontsize=20)
        #plt.ylabel('Flux')
        plt.suptitle('Kepler ' + keplerID,fontsize=20,fontweight='bold')
        #plt.figlegend((axs[0],axs[1]),('GALEX','Kepler'),'upper left')
        f.text(0.1, 0.5, 'Flux', ha='center', va='center', rotation='vertical',fontsize=20)

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        
        plt.show()
