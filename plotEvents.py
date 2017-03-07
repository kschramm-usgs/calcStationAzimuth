#!/bin/env python

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

client = Client("IRIS")

starttime = UTCDateTime("2017-02-01")
endtime = UTCDateTime("2017-03-07")

cat = client.get_events(starttime=starttime, endtime=endtime,
                        minmagnitude=5)
print(cat)

cat.plot()
