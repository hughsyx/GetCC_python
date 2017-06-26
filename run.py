from GetANIR import *
import time
from obspy.core import Trace, Stream, AttribDict
jt = time.time()
sta1 = 'CI.PSD'
sta2 = 'CI.SNCC'
comp1 = 'BHZ'
comp2 = 'BHZ'
ANIR = MyGetANIR(sta1,sta2,comp1,comp2)
print ">>> Job Finished. It took %.2f seconds" %(time.time()-jt)
plt.plot(ANIR)
plt.show()

pair = "%s:%s" % (sta1,sta2)
out_path = '/Users/Hugh/Documents/C1_results/'+sta1+'/'
filename = out_path + sta1.split('.')[1] + '_' + sta2.split('.')[1]+'.SAC'
try:
	os.makedirs(os.path.split(filename)[0])
except:
	pass

mytrace = Trace(data=ANIR)
mytrace.stats['station'] = pair
mytrace.stats.sac = AttribDict()
mytrace.stats.sac.npts = len(ANIR)

st = Stream(traces=[mytrace,])
st.write(filename,format ='SAC')
del st
