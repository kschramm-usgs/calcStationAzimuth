#!/bin/env python

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics.base import locations2degrees
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel
from numpy import sin, cos

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
inventory = client.get_stations(network="IU", station="ANMO", 
                                channel="BH1",level="response",
                                starttime=eventTime)

# next, get the station coordinates
print("getting station coordinates")
station_coordinates = []
for network in inventory:
    for station in network:
        for channel in station:
            station_coordinates.append((network.code, station.code, 
                                        station.latitude, station.longitude, 
                                        station.elevation,channel.azimuth))

# then for each station in the list get the distance and azimuth
# need to think about what source-receiver distances we want to use
# for p-waves.
# pick a model for estimating arrival times
print("calculating travel times and requesting data")
model = TauPyModel(model="iasp91")
for station in station_coordinates:
# first calculate the source-receiver distance
    DegDist = locations2degrees(eventLat, eventLon,
                                station[2], station[3])
    StationAziExpec = gps2dist_azimuth(eventLat, eventLon,
                                       station[2], station[3]) 
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
# ask for data 200 s before and 50 s after P time
# the larger window at the beginning is so that we can look at the SNR
        bTime=arrTime-200
        eTime=arrTime+50
        try:
            st = client.get_waveforms(station[0],station[1],"00","BH?",
                                      bTime,eTime,attach_response=True)
        except:
            print("No data for station "+station[1])
	    continue #use a continue to go back to the beginning of the loop
#        prefilt = (1/4.,1/2., 10., 20.) # this may need to be changed, but 
#                                         I will need to discuss with tyler
#        st.remove_response(output="DISP",pre_filt=prefilt)
# use this line if you want to see the plots
#        st.remove_response(output="DISP",pre_filt=prefilt,plot=True)
# Break up the stream into traces to remove the gain
        BH1 = st[0]
        BH2 = st[1]
        BHZ = st[2]
        print('check the component')
        print('check the component')
        BHZ.plot()
        BH1.plot()

        st[0] = BH1.remove_sensitivity()
        st[1] = BH2.remove_sensitivity()
        st[2] = BHZ.remove_sensitivity()
         
# take a look at the data
        st.plot()
# make sure we are in NE orientation
        BHN = st[0]
        BHE = st[1]
        BHZ = st[2]
        station_orientation = station[5]
        if (station_orientation != 0.0):
            BHE.data = ( sin(station_orientation)*st[0] 
                +   cos(station_orientation)*st[1])
            BHN.data = (-sin(station_orientation)*st[2] 
                +   cos(station_orientation)*st[0])
        BHE.plot()
        BHN.plot()

        
# next we need to filter.
#        tr_filt = BHZ
#        tr_filt.filter('lowpass', freq=0.5, corners = 4, zerophase = True) 
         
         
    else:
        print("Station "+ station[1] +"doesn't fit in parameters for P-wave arrivals")

