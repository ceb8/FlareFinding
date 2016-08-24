#!/usr/bin/env python

from __future__ import print_function

from re import search

from astropy.coordinates import SkyCoord


def makeObjectString(objLst):
    outLst = [objLst[x] for x in [0,5,2,3]]
    
    coord = SkyCoord(b = float(objLst[6]),l = float(objLst[7]),
                     unit='deg', frame='galactic')
    outLst.append(str(coord.icrs.ra.deg))
    outLst.append(str(coord.icrs.dec.deg))
    
    return ','.join(outLst)
              


# Open list of Kepler observed objects
KEPLER = open('apjs489745t4_mrt.txt','r')

# data does not start until line 34
for _ in xrange(34):
    next(KEPLER)

#Extract Kepler IDs and put in a list
keplerIdList = []
for obs in KEPLER:
    keplerIdList.append(search('([0-9]{6,8}) ', obs).group(1))

print("Number of Kepler targets: ",len(keplerIdList))

# we are done with the Kepler file
KEPLER.close()

# Open list of GALEX objects in the Kepler Input Catalog (csv)
KGMATCH = open('KGMatch.csv','r')

# first line just tells what each column is
next(KGMATCH)

# Look for each GALEX entry in the Kepler list
# and if found record relevant data
KGOverlap = ["kic_kepler_id,glx_objid,kic_deg_ra,kic_deg_dec,glx_deg_ra,glx_deg_dec"]
i = 0
for obs in KGMATCH:
    obsLst = obs.split(',')
    if obsLst[0] in keplerIdList:
        KGOverlap.append(makeObjectString(obsLst))
            
    if not i%1000:
        print(i)
    i += 1

print("Number of KGMatch entries: ",i)
    
# we are done with KGMatch
KGMATCH.close()

# outputing our results
print('The number of kepler targets also targeted by GALEX is: ',
      len(KGOverlap)-1)

KGOUT = open('KGCommonObservations.csv','w')
for line in KGOverlap:
    KGOUT.write(line + '\n')
KGOUT.close()
