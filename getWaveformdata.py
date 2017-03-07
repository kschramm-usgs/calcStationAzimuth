#!/bin/env/ python

''' 
Routine to get waveform data from IRIS given network and station
and event information
'''

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.taup import TauPyModel

# set the client where you want to retrieve the data
client = Client("IRIS")

# event time - this can be calculate based on the p-wave arrival
# this is for an event in Bolivia
t = UTCDateTime("2017-02-21T14:09:04.000")
evMag = 6.5
evDepth = 597.9
evLat = -19.284
evLon = 63.899

stLat = 34.945910
stLon = -106.457200
stDepth = -1820

arrivals = model.get_travel_times_geo(source_depth_in_km=evDepth,
                                      source_latitude_in_deg = evLat,
                                      source_longitude_in_deg = evLon,
                                      receiver_latitude_in_deg = stLat,
                                      receiver_longitude_in_deg = stLon)
