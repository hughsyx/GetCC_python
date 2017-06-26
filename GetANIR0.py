import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
from obspy.signal.invsim import waterlevel
import pdb
#import traceback
#import glob
import time
#import calendar
import sys
from obspy.core import Trace, Stream, AttribDict
from obspy.core import utcdatetime, UTCDateTime,read
import os
from obspy.core import utcdatetime, UTCDateTime,read
from CANIR import *
import gc
""" this function is used to get the ambient noise impulse response """

def myReadData(sta1,sta2,comp1,comp2,day):
	
	data_path1 = '/Users/Hugh/Documents/One\ Drive/OneDrive/data/%s/%s.*.%s.%s.SAC' % (comp1,sta1.replace('.','_'),sta1,comp1,day)
        data_path2 = '/Users/Hugh/Documents/One\ Drive/OneDrive/data/%s/%s.*.%s.%s.SAC' % (comp2,sta2.replace('.','_'),sta2,comp2,day)
	#pdb.set_trace()
	st1 = read(data_path1,dytpe=np.float)
	st2 = read(data_path2,dytpe=np.float)
	del data_path1, data_path2
	return st1,st2	

def CheckFile(sta1,sta2,comp1,comp2,day):
	file_exist = 0
    #pdb.set_trace()
	data_path1 = '/Users/Hugh/Documents/One\ Drive/OneDrive/data/%s/%s..%s.%s.SAC' % (comp1,sta1.replace('.','_'),sta1,comp1,day)
	data_path12 = '/Users/Hugh/Documents/One\ Drive/OneDrive/data/%s/%s.00.%s.%s.SAC' % (comp1,sta1.replace('.','_'),sta1,comp1,day)
        data_path2 = '/Users/Hugh/Documents/One\ Drive/OneDrive/data/%s/%s..%s.%s.SAC' % (comp2,sta2.replace('.','_'),sta2,comp2,day)
	data_path22 = '/Users/Hugh/Documents/One\ Drive/OneDrive/data/%s/%s.00.%s.%s.SAC' % (comp2,sta2.replace('.','_'),sta2,comp2,day)
        data_path23 = '/Users/Hugh/Documents/One\ Drive/OneDrive/data/%s/%s.10.%s.%s.SAC' % (comp2,sta2.replace('.','_'),sta2,comp2,day)
        data_path13 = '/Users/Hugh/Documents/One\ Drive/OneDrive/data/%s/%s.10.%s.%s.SAC' % (comp1,sta1.replace('.','_'),sta1,comp1,day)
	if (os.path.isfile(data_path1)+os.path.isfile(data_path12)+os.path.isfile(data_path13)) & (os.path.isfile(data_path2)+os.path.isfile(data_path22)+os.path.isfile(data_path23)):
		file_exist = 1
	del data_path1,data_path2
	return file_exist	
	
def Day2String(day):
	if day < 10:
		Day = '00'+str(day)
	elif day < 100:
		Day = '0'+str(day)
	else:
		Day = str(day)	
	del day
	return Day

def MyGetANIR(sta1,sta2,comp1,comp2,outpath): 

	dt = 0.2 # 1/sampling rate
	winlen = int(30*60/dt) # in seconds
	winover = int(20*60/dt) # in seconds
	#nmb_w = 0
	i = 0
	eps = 0.001
	npt = 10
	cc = 0
	N = 0
	ANIR = 0
	outpath2 = outpath + '/' + comp1 + '_' + comp2 + '/' + sta1 + '/' + sta2 + '/'
	for i in range(0,365):
		day = i+1
		Day = Day2String(day)
		#pdb.set_trace()
		if CheckFile(sta1,sta2,comp1,comp2,Day)==0 :
			continue	
		st1,st2 = myReadData(sta1,sta2,comp1,comp2,Day)
		tr1 = st1[0]
		tr2 = st2[0]
		v1 = tr1.filter("bandpass",freqmax = 1, freqmin = 0.1).data
		v2 = tr2.filter("bandpass",freqmax = 1, freqmin = 0.1).data
		max1 = max(abs(v1))
		max2 = max(abs(v2))
		Len = min(len(v1),len(v2))
		nwin = int(np.floor((Len-winlen)/winover)+1)
		nmb_w = 0
		temp = 0
		filename = outpath2 + Day + '.SAC'
		del st1, st2
		for j in range(nwin):
			dur = range(j*winover,j*winover+winlen)
			v12 = v1[dur]
			v22 = v2[dur]
			if (max(abs(v12)) < np.std(v12)*4) * (max(abs(v22)) < np.std(v22)*4):
				if (sum(abs(v12)<eps*max1) < 500) * (sum(abs(v22)<eps*max2) < 500):
					nmb_w = nmb_w + 1
					c=myCorr2(v22,v12,npt)
					temp = temp + c
		if nmb_w >20:
			cc = cc+temp
			N = N+nmb_w
		gc.collect()
	if N > 0:
		ANIR = cc/N
	return ANIR,N 
		
def export_sac(c,filename,sta1,sta2,nmb_w):
   	try:
        	os.makedirs(os.path.split(filename)[0])
    	except:
        	pass
    	mytrace = Trace(data=c)
    	mytrace.stats['station'] = sta1 + '_' + sta2
    	mytrace.stats.sac = AttribDict()
    	mytrace.stats.sac.npts = len(c)
	mytrace.stats.starttime = nmb_w
    	st = Stream(traces=[mytrace,])
    	st.write(filename,format = 'SAC')
    	del st




















