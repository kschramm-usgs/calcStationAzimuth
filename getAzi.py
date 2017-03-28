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
import sys
import argparse
import os


''' 
This function will calculate the azimuth of a station from particle motion.

This is a translation of Adam Ringler's matlab function.  

'''

def getargs():
    """ Grab command line arguments to run program. """
    parser = argparse.ArgumentParser(description = "Program to calculate station azimuth from p-wave particle motion")

    parser.add_argument('-resDir',type=str, action = "store",
                        dest="resDir", required=True,
                        help="Result directory name Example: blah")

    parser.add_argument('-eventTime', type=str, action = "store",
                         dest="eventTime", required=True,
                         default ="2017-02-21T14:09:04",
                         help="Enter a UTC time, ie: 2017-02-21T14:09:04")

    parser.add_argument('-eventLat',type=float, action = "store",
                         dest="eventLat", required=True,
                         default =-19.284,
                         help="Enter the lat, lon, depth for the event")

    parser.add_argument('-eventLon',type=float, action = "store",
                         dest="eventLon", required=True,
                         default =-63.899,
                         help="Enter the longitude for the event")

    parser.add_argument('-n', type=str, action="store",
                         dest = "network", required=True,
                         help="Network name Example: IU")

    parser.add_argument('-sta', type=str, action="store",
                        dest="sta", required = False,
                        help="Stations to use. Example with a \
                                comma (,) separator : TUC,ANMO")

    parser.add_argument('-cha', type=str, action="store",
                        dest="cha", required = False,
                        help="Channels to use. Example: BH*")

    parserval=parser.parse_args()
    return parserval

# Start of the main part of the program
if __name__ == "__main__":

#use a modified version of AdamR's parser function to get information 
    parserval = getargs()
# event in bolivia
    eventTime = UTCDateTime(parserval.eventTime)
    eventLat = parserval.eventLat
    eventLon = parserval.eventLon
    eventDepth = 597

#check for existence of the result directory
    resultdir = parserval.resDir
    if resultdir[-1] == '/':
        resultdir = resultdir[:-1]

    if not os.path.exists(os.getcwd() + '/' + resultdir):
        os.mkdir(os.getcwd() + '/' + resultdir)

    statfile = open(os.getcwd()+'/'+parserval.resDir+'/Results_'+
                    str(eventTime)+'_'+parserval.network+'.csv','w')
    statfile.write('station, channel, expected Baz, calc Baz 1, calcBaz2, linearity\n')

# event in philipines
#eventTime = UTCDateTime("2017-01-10T06:13:48")
#eventLat = 4.478
#eventLon = 122.617
#eventDepth = 627.2

    net = parserval.network
    stat = parserval.sta
    chan = parserval.cha

# get station info
# first build the inventory of stations
    print("building station list for the event")
    client = Client("IRIS")
    inventory = client.get_stations(network=net, station=stat,
                                    channel=chan,level="response",
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
    model = TauPyModel(model="ak135")
    for station in station_coordinates:
# first calculate the source-receiver distance
        DegDist = locations2degrees(eventLat, eventLon,
                                    station[2], station[3])
# need to add tolerance for distance so that we are only using P-arrivals
# need to talk to tyler about which P should be used?  P? Pdiff? pP? PP?
# tyler also feels we really want a direct P at teleseismic distances, so 
# let us start with 25-90
        if DegDist > 25 and DegDist < 90:
            print("Station "+station[1]+" will have a P-wave arrival")
            StationAziExpec = gps2dist_azimuth(eventLat, eventLon,
                                               station[2], station[3]) 
            print("station lat, lon:"+str(station[2])+","+str(station[3]))
            statBaz = StationAziExpec[2]
            print "The expected back azimuth for "+station[1]+" is "+str(statBaz)
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
            #st.plot()

# make sure we are in NE orientation
# we may decide to take this step out - find out if we assume N/S does
# the program spit out the azimuth we have in the metadata?
            BHN = st[0].copy()
            BHE = st[1].copy()
            BHZ = st[2].copy()
            station_orientation = station[5]
            if (station_orientation != 0.0):
                statO = np.radians(station_orientation)
                print("Rotating station to N/E. Correction Angle is: "
                      +str(station_orientation))
                #BHN.data = (sin(station_orientation)*BHN.data
                #    +   cos(station_orientation)*BHE.data)
                #BHE.data = (sin(-station_orientation)*BHE.data
                #    +   cos(station_orientation)*BHN.data)
                BHN.data = (sin(statO)*BHN.data
                    +   cos(statO)*BHE.data)
                BHE.data = (sin(-statO)*BHE.data
                    +   cos(statO)*BHN.data)
                
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
            SignalStart = arrTime-1
            SignalEnd = arrTime+8
            SignalBHN.trim(SignalStart, SignalEnd)
            SignalBHE.trim(SignalStart, SignalEnd)
            SignalBHZ.trim(SignalStart, SignalEnd)
     
            SNR_BHN = (SignalBHN.std()**2)/(NoiseBHN.std()**2)
            SNR_BHE = (SignalBHE.std()**2)/(NoiseBHE.std()**2)
            SNR_BHZ = (SignalBHZ.std()**2)/(NoiseBHZ.std()**2)
     
            #print(SNR_BHN, SNR_BHE, SNR_BHZ)
# normalize
            BHNmax=np.max(abs(SignalBHN.data))
            BHNmin=np.min(abs(SignalBHN.data))
            BHEmax=np.max(abs(SignalBHE.data))
            BHEmin=np.min(abs(SignalBHE.data))
            normFac=np.max([BHNmax,BHNmin,BHEmax,BHEmin])
            SignalBHN.data = SignalBHN.data/normFac
            SignalBHE.data = SignalBHE.data/normFac
            SignalBHZ.data = SignalBHZ.data/normFac
            plt.subplot(3,1,1)
            plt.suptitle('Filtered, rotated and cut waveforms for %s'%(station[1]))
            plt.plot(SignalBHN.data,label='BHN')
            plt.legend()
            #SignalBHN.plot()
            plt.subplot(3,1,2)
            plt.plot(SignalBHE.data,label='BHE')
            plt.legend()
            plt.subplot(3,1,3)
            plt.plot(SignalBHZ.data,label='BHZ')
            plt.legend()
            plt.show()
    
# time to get serious!  we are ready to do the actual calculation!!!!!!!!
            A = np.transpose(np.matrix(SignalBHE.data))
            b = SignalBHN.data
            #print A.shape,b.shape

            lresult = lsqr(A,b)
            ang = np.degrees(np.arctan2(1.,lresult[0]))
            #print "The answer you are looking for:"
            #print(ang)

# Adam uses this to calculate the linearity.
            BHNsq = sum(SignalBHN.data*SignalBHN.data)
            BHNEsq = sum(SignalBHN.data*SignalBHE.data)
            BHEsq = sum(SignalBHE.data*SignalBHE.data)
            #print BHNsq, BHNEsq, BHEsq
            eigMat = np.matrix([[BHNsq, BHNEsq], [BHNEsq, BHEsq]])
            eigd,eigv = eig(eigMat)
            #print "The eigen vals:"
            #print eigd
            #print "The eigen vecs:"
            #print eigv.real
            line = (eigd[1]/(eigd[0]+eigd[1]))-(eigd[0]/(eigd[0]+eigd[1]))       
            #print "The linearity:"
            #print line

# what is the motivation for these calculations?
# now do some stuff about the quadrant?
            ang2 = np.degrees(np.arctan2(eigv[0][1],eigv[1][1]))
            #print "This is a different value for the angle:"
            #print ang2
            #print ang
            if abs(statBaz-(ang2+180))<abs(statBaz-ang2):
                ang2 = np.mod(ang2+180,360)
                print"ang2 is 180 off: "+str(ang2) 
            if (ang2<0):
                ang2 = ang2+180
                print"ang2 lt 0"
            if(abs(statBaz-(ang+180))<abs(statBaz-ang)):
                ang=np.mod(ang+180,360)
                print"ang is 180 off: "+str(ang) 
            if(ang<0):
                ang=ang+180;
                print"ang lt 0"
            print("The calculated values are: "+str(ang)+" and "+str(ang2))
# capture some statistics
# there is likely a better way to do this, but in the interest of time...
#            statfile.write('%s,%s,%s,%s,%s\n'%
#                    station[1],station[2],statBaz,ang,ang2,line)
            statfile.write(station[0]+', ')
            statfile.write(station[1]+', ')
            statfile.write(str(statBaz)+', ')
            statfile.write(str(ang)+', ')
            statfile.write(str(ang2)+', ')
            statfile.write(str(line)+'\n')

# don't look at the plot if it's not worth your time.
            if abs(line) < 0.8:
                print("Linearity value bad.  Skipping station.")
                continue
# now create a nice plot. 
            ax = plt.subplot(111, projection='polar')
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            theta = np.arctan2(SignalBHE.data,SignalBHN.data)
            r = np.sqrt(SignalBHE.data*SignalBHE.data 
                  + SignalBHN.data*SignalBHN.data)
            calcR = [1.5, 1.5]
            calcTheta = [np.radians(ang),np.radians(ang+180)]
            calcTheta2 = [np.radians(ang2),np.radians(ang2+180)]
            expcTheta = ([np.radians(StationAziExpec[2]), 
                          np.radians(StationAziExpec[2]+180)])
            label1='Calculated Baz = '+str(ang)
            label2='Calculated Baz 2 = '+str(ang2)
            label3='Expected Baz = '+str(StationAziExpec[2])
            plt.plot(calcTheta,calcR,'blue',label=label1)
            plt.plot(calcTheta2,calcR,'cyan',label=label2)
            plt.plot(expcTheta,calcR,'black',label=label3)
            plt.plot(theta,r,'red',label='Particle Motion')
            plt.legend(bbox_to_anchor=(0., 1.02, 1., 0.102),loc=3,borderaxespad=0.)
            plt.text(np.pi,2,str(station[1]))
            plt.show()



        else:
            print("Station "+ station[1] +" doesn't fit in parameters for P-wave arrivals")
