#!/bin/env python

from obspy.geodetics.base import locations2degrees
# for vmqc001 you need the import below.
#from obspy.core.util import locations2degrees

staLat = 34.945910 # default station is ANMO
staLon = -106.4572 

eqLat = 37.580 #default event in Turkey, 2017055 11:07:27, M5.6
eqLon = 38.440

# get the distance between the earthquake and station in degrees
DegDist = locations2degrees(staLat, staLon, eqLat, eqLon)

# print the result
print 'The distance between the station and earthquake is: '+str(DegDist)
