#!/bin/env python

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics.base import locations2degrees
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel
from numpy import sin, cos
from numpy import arctan as atan
from scipy.sparse.linalg import lsqr
from scipy.linalg import eig
import numpy as np
import matplotlib.pyplot as plt


''' 
This function will calculate the azimuth of a station from particle motion.

This is a translation of Adam Ringler's matlab function.  

'''

# event in bolivia
#eventTime = UTCDateTime("2017-02-21T14:09:04")
#eventLat = -19.284
#eventLon = -63.899
#eventDepth = 597

# event in philipines
eventTime = UTCDateTime("2017-01-10T06:13:48")
eventLat = 4.478
eventLon = 122.617
eventDepth = 627.2

#station = "ANMO"
#network = "IU"

# get station info
# first build the inventory of stations
print("building station list for the event")
client = Client("IRIS")
inventory = client.get_stations(network="IU", station="CTAO", 
                                channel="BH1",level="response",
                                location="00",starttime=eventTime)

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
    print "The expected back azimuth:"
    print(StationAziExpec)
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

        st[0] = BH1.remove_sensitivity()
        st[1] = BH2.remove_sensitivity()
        st[2] = BHZ.remove_sensitivity()
         
# take a look at the data
        st.plot()

# make sure we are in NE orientation
# we may decide to take this step out - find out if we assume N/S does
# the program spit out the azimuth we have in the metadata?
        BHN = st[0]
        BHE = st[1]
        BHZ = st[2]
        station_orientation = station[5]
        if (station_orientation != 0.0):
            BHE.data = ( sin(station_orientation)*st[0] 
                +   cos(station_orientation)*st[1])
            BHN.data = (-sin(station_orientation)*st[2] 
                +   cos(station_orientation)*st[0])
        
# next we need to filter.
        BHN.detrend('demean')
        BHE.detrend('demean')
        BHZ.detrend('demean')
        
        BHN.taper(max_percentage=0.05)
        BHE.taper(max_percentage=0.05)
        BHZ.taper(max_percentage=0.05)
# now, we actually filter!
        BHN.filter('lowpass', freq=0.5, corners=4, zerophase=True)
        BHE.filter('lowpass', freq=0.5, corners=4, zerophase=True)
        BHZ.filter('lowpass', freq=0.5, corners=4, zerophase=True)
         
# calculate the SNR - I forsee putting in something here to skip calculation
# if the SNR is below a defined tolerance.
        NoiseBHN = BHN.copy()
        NoiseBHE = BHE.copy()
        NoiseBHZ = BHZ.copy()
        NoiseStart = BHN.stats.starttime + 10
        NoiseEnd = BHN.stats.starttime + 150
        NoiseBHN.trim(NoiseStart, NoiseEnd)
        NoiseBHE.trim(NoiseStart, NoiseEnd)
        NoiseBHZ.trim(NoiseStart, NoiseEnd)
     
# need to try other times than iasp.  austin says usgs uses ak135
        SignalBHN = BHN.copy()
        SignalBHE = BHE.copy()
        SignalBHZ = BHZ.copy()
        SignalStart = arrTime
        SignalEnd = arrTime+10
        SignalBHN.trim(SignalStart, SignalEnd)
        SignalBHE.trim(SignalStart, SignalEnd)
        SignalBHZ.trim(SignalStart, SignalEnd)
     
        SNR_BHN = (SignalBHN.std()**2)/(NoiseBHN.std()**2)
        SNR_BHE = (SignalBHE.std()**2)/(NoiseBHE.std()**2)
        SNR_BHZ = (SignalBHZ.std()**2)/(NoiseBHZ.std()**2)
     
        print(SNR_BHN, SNR_BHE, SNR_BHZ)
# normalize
# the obspy does not seem to be working...
#        BHN.normalize
#        BHE.normalize
#        BHZ.normalize
#        BHN.plot()
#        BHE.plot()
#        BHZ.plot()
#normalize with Adam's method
#        BHNmax=SignalBHN.max()
#        BHNmax=np.max(SignalBHN.data)
#        print "bhn max"
#        print (BHNmax, BHNmaax)
#        BHNmin=np.min(SignalBHN.data)
#        print (BHNmin)
#        BHEmax=np.max(SignalBHE.data)
#        BHEmin=np.min(SignalBHE.data)
#        BHN_range=BHNmax-BHNmin
#        BHE_range=BHEmax-BHEmin
#        BHHmax=np.max([BHN_range,BHE_range])
#        BHNmine=BHHmax*BHNmin/BHN_range
#        BHEmine=BHHmax*BHEmin/BHE_range
#        BHNmaxe=BHHmax*BHNmax/BHN_range
#        BHEmaxe=BHHmax*BHEmax/BHE_range
#        SignalBHN.data=SignalBHN.data+(BHNmaxe+BHNmine+BHNmax-BHNmin)/2
#        SignalBHE.data=SignalBHE.data+(BHEmaxe+BHEmine+BHEmax-BHEmin)/2
#        SignalBHN.plot()
#        SignalBHE.plot()
# ok, that also didn't work.  grrrr. why is this so complicated.  now
# we try the Kim Schramm method.    
        BHNmax=np.max(abs(SignalBHN.data))
        BHNmin=np.min(abs(SignalBHN.data))
        BHEmax=np.max(abs(SignalBHE.data))
        BHEmin=np.min(abs(SignalBHE.data))
        SignalBHN.data = SignalBHN.data/BHNmax
        SignalBHN.plot()
        SignalBHE.data = SignalBHE.data/BHEmax
        SignalBHE.plot()
        


    
# time to get serious!  we are ready to do the actual calculation!!!!!!!!
        A = np.transpose(np.matrix(SignalBHE.data))
        b = SignalBHN.data
        print A.shape,b.shape

        lresult = lsqr(A,b)
        ang = np.degrees(np.arctan2(1.,lresult[0]))
        print "The answer you are looking for:"
        print(ang)

# Adam uses this to calculate the linearity.
        BHNsq = sum(SignalBHN.data*SignalBHN.data)
        BHNEsq = sum(SignalBHN.data*SignalBHE.data)
        BHEsq = sum(SignalBHE.data*SignalBHE.data)
        print BHNsq, BHNEsq, BHEsq
        eigMat = np.matrix([[BHNsq, BHNEsq], [BHNEsq, BHEsq]])
        eigd,eigv = eig(eigMat)
        print "The eigen vals:"
        print eigd
        print "The eigen vecs:"
        print eigv
        line = (eigd[1]/(eigd[0]+eigd[1]))-(eigd[0]/(eigd[0]+eigd[1]))       
        print "The linearity:"
        print line

# what is the motivation for these calculations?
# now do some stuff about the quadrant?
        ang2 = np.degrees(np.arctan2(eigv[0][1],eigv[1][1]))
        print "This is a different value for the angle:"
        print ang2
        if (ang2<0):
            ang2 = ang2+180
        if abs(StationAziExpec[1]-(ang2+180))<abs(StationAziExpec[1]-ang2):
            ang2 = np.mod(ang2+180,360)
        if(ang<0):
            ang=ang+180;
        if(abs(StationAziExpec[1]-(ang+180))<abs(StationAziExpec[1]-ang)):
            ang=np.mod(ang+180,360)
        print "This is a different angle:"
        print ang2
        print "Not sure why Adam is calculating this:"
        print(ang)
# what about plotting particle motion...
# need to get amplitudes in plot them up.
# would be nice to see the a line at the calculated angle...
# can we look at the phase of the seismograms to get the quadrant?
# now create a nice plot. 
        ax = plt.subplot(111, projection='polar')
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        theta = np.arctan2(SignalBHE.data,SignalBHN.data)
        r = np.sqrt(SignalBHE.data*SignalBHE.data 
                  + SignalBHN.data*SignalBHN.data)
        #print r
        #print theta
        plt.plot(theta,r,'red',label='Particle Motion')
        calcR = [1.5, 1.5]
        calcTheta = [np.radians(ang),np.radians(ang+180)]
        expcTheta = ([np.radians(StationAziExpec[2]), 
                      np.radians(StationAziExpec[2]+180)])
        plt.plot(calcTheta,calcR,'blue',label='Calculated Baz')
        plt.plot(expcTheta,calcR,'black',label='Expected Baz')
        ax.legend(loc='upper left')
        plt.show()


    else:
        print("Station "+ station[1] +"doesn't fit in parameters for P-wave arrivals")

