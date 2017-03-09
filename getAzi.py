#!/bin/env python

from obspy import UTCDateTime

''' 
This function will calculate the azimuth of a station from particle motion.

This is a translation of Adam Ringler's matlab function.  

'''

eventTime = UTCDateTime("2017-02-21T14:09:04")
eventLat = -19.284
eventLon = -63.899

station = "ANMO"
network = "IU"

# get station info
