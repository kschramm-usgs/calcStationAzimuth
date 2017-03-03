#!/bin/env python
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

eventLat = -19.284
eventLon = -63.899
starttime = UTCDateTime("2017-02-21")
endtime = UTCDateTime("2017-02-22")

client = Client("IRIS")

inventory = client.get_stations(network ="IU,US,CU,IW,II", latitude=eventLat, longitude = eventLon, maxradius = 20, starttime = starttime, endtime = endtime)

print(inventory)
