# calcStationAzimuth
This project will calculate the azimuth between source receiver pairs.

It uses similar command line input to Adam Ringler's syncomp.py.  Here is an example:

./getAzi.py -resDir junk -n IU -cha "BH1" -eventTime '2017-02-21T14:09:04' -eventLat -19.284 -eventLon -63.899 -eventDepth 597 -sta "*"

This is for an event in Bolivia.

Need to have commandLineInfo from waveformUtils repository to run.
