#!/usr/bin/env python

from __future__ import print_function

from gPhoton.gFind import gFind
from astropy.coordinates import SkyCoord

from time import gmtime

from datetime import date
from time import mktime

###########################
#
# Functions
#
###########################


# Test whether a time given as GALEX time is within
# the Kepler observation time range (2009-2013)
# 0 == within Kepler range
# -1 == before Kepler commenced observations 
# 1 = after Kepler ceased observations
def compareKeplerRange(galexTime):
    
    # "GALEX Time" = "UNIX Time" - 315964800
    utcTime = gmtime(galexTime + 315964800)
    
    #print(utcTime.tm_year)

    if utcTime.tm_year < 2009:
        return -1
        
    if utcTime.tm_year > 2013:
        return 1

    return 0
    

# Takes a gFind result and returns the observation times
# as a string in the following format (NUV,FUV):
# ('(t0_0:t0_1),(t1_0,t1_1),...','(t0_0:t0_1),(t1_0,t1_1),...')
def gObsTimes2Str(galexObs, nuv = True, fuv = False):
    nuvTimes = ""
    fuvTimes = ""
    
    if nuv:
        if 'NUV' in galexObs:
            initTimes = galexObs['NUV']['t0']
            finTimes = galexObs['NUV']['t1']
            times = zip(initTimes,finTimes)
            nuvTimes = str(times)[1:-1].replace("), ","),").replace(", ",":")
        
    if fuv:
        if 'FUV' in galexObs:
            initTimes = galexObs['FUV']['t0']
            finTimes = galexObs['FUV']['t1']
            times = zip(initTimes,finTimes)
            nuvTimes = str(times)[1:-1].replace("), ","),").replace(", ",":")

    return (nuvTimes,fuvTimes)

###########################
#
# Main program
#
###########################

# Go through each of the extracted Kepler/GALEX overlap
# Check the total observation time is > 30min
# if yes extract obs datetimes and check if it's in Kepler date range
# (this mean figuring out how to get a real datetime out of galex datetimes)
# if yes output to new file
# info in output file:
# everything in input file, plus GALEX time info


# Open csv containing objects observed by Kepler and Galex
KGOBJ = open('KGCommonObservations.csv','r')

# The first line just tells what each column is
next(KGOBJ)
# with added appending, skip to last done entry
for _ in xrange(42807):
    next(KGOBJ)

# Open outfile (appending)
KGOUT = open('KGSimultaneousObservations.csv','a')

if not KGOUT.tell(): # if file is empty, write header
    KGOUT.write("kic_kepler_id,glx_objid,kic_deg_ra,kic_deg_dec,glx_deg_ra,glx_deg_dec,(G_t0:G_t1),...\n")


i = 0
j = 0
k = 42807
l = 27664
for obj in (KGOBJ):

    if k%50 == 0:
        print(k)
    k+=1

    obj = obj.strip('\n')
    objLst = obj.split(',')

    # turn the given lat long coordiante into ra and dec
    #objCoord = SkyCoord(b = float(objLst[5]),l = float(objLst[6]),
    #                     unit='deg', frame='galactic')
    #icrsCoord = objCoord.icrs
    
    # get GALEX obsevation data for this object
    galexObs = gFind(band='NUV',
                     #skypos=[icrsCoord.ra.deg,icrsCoord.dec.deg],
                     skypos=[float(objLst[4]),float(objLst[5])],
                     minexp=300,
                     trange=[914821200,1072501200], # 2009 - 20013, kepler range
                     quiet=True)

    # check total time
    totalTime = galexObs['NUV']['expt']
    #print(totalTime)
    if totalTime < 1800: # observation time under 30 min, go to next entry
        i+=1
        continue
  
    # check if the observation time is in Keplers range
    # Fails if the end of the last observation is < start of Kepler obs
    # or start of first observation is > end of Kepler obs
    
    #if (( compareKeplerRange(galexObs['NUV']['t1'][-1]) == -1 ) or
    #    ( compareKeplerRange(galexObs['NUV']['t0'][0]) == 1)):
    #    j+=1
    #    continue

    # At this point we have a simultaneous Kepler/Galex obsevation
    # So save off the information
    l+=1
    KGOUT.write(obj + "," + gObsTimes2Str(galexObs)[0] + '\n')




KGOBJ.close()
    
# outputing our results
print("The number of GALEX observations that were too short is: ", i)
#print("The number of longer GALEX observation not simultaneous to Kepler is: ",j)
print('The number ofsumultaneos Kepler and GALEX targets is: ', l)


#close the outfile
KGOUT.close()


    
    
    

    
    
    
