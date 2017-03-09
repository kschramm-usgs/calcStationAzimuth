#!/bin/env python

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics.base import locations2degrees
from obspy.taup import TauPyModel

''' 
This function will calculate the azimuth of a station from particle motion.

This is a translation of Adam Ringler's matlab function.  

'''

eventTime = UTCDateTime("2017-02-21T14:09:04")
eventLat = -19.284
eventLon = -63.899
eventDepth = 597

station = "ANMO"
network = "IU"

# get station info
# first build the inventory of stations
inventory = client.get_stations(network="IU", station="*", 
                                starttime=eventTime)

# next, get the station coordinates
station_coordinates = []
for network in inventory:
    for station in network:
        station_coordinates.append((network.code, station.code, 
                                    station.latitude, station.longitude, 
                                    station.elevation))

# then for each station in the list get the distance and azimuth
# need to think about what source-receiver distances we want to use
# for p-waves.
for station in station_coordinates:
    DegDist = locations2degrees(station[3], station[4], eventLat, eventLon)
# need to add tolerance for distance so that we are only using P-arrivals
# need to talk to tyler about which P should be used?  P? Pdiff? pP? PP?
    arrivals = model.get_travel_times(source_depth_in_km = eventDepth
                                      distance_in_degree=DegDist)
    

