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
print("building station list for the event")
client = Client("IRIS")
inventory = client.get_stations(network="IU", station="*", 
                                starttime=eventTime)

# next, get the station coordinates
print("getting station coordinates")
station_coordinates = []
for network in inventory:
    for station in network:
        station_coordinates.append((network.code, station.code, 
                                    station.latitude, station.longitude, 
                                    station.elevation))

# then for each station in the list get the distance and azimuth
# need to think about what source-receiver distances we want to use
# for p-waves.
# pick a model for estimating arrival times
print("calculating travel times and requesting data")
model = TauPyModel(model="iasp91".)
for station in station_coordinates:
    DegDist = locations2degrees(station[2], station[3], eventLat, eventLon)
# need to add tolerance for distance so that we are only using P-arrivals
# need to talk to tyler about which P should be used?  P? Pdiff? pP? PP?
# tyler also feels we really want a direct P at teleseismic distances, so 
# let us start with 25-90
    if DegDist > 25 and DegDist < 90:
        arrivals = model.get_travel_times(source_depth_in_km = eventDepth,
                                          distance_in_degree=DegDist,
                                          phase_list = ["P"])
# now use the arrival time to get some data
        arrTime=eventTime + arrivals[0].time
# ask for data one minute before and 5 minutes after P time
        bTime=arrTime-60
        eTime=arrTime+300
        try:
            st = client.get_waveforms(station[0],station[1],"00","BHZ",
                                      bTime,eTime)
        except:
            print("No data for station "+station[1])

        st.plot()

    else:
        print("Station "+ station[1] +"doesn't fit in parameters for P-wave arrivals")

    

