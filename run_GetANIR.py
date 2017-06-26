from GetANIR import *
import time
from obspy.core import Trace, Stream, AttribDict
from pyrocko import model
import pdb
import warnings
warnings.filterwarnings("ignore")

def run(sta1,sta2,comp1,comp2):
	jt = time.time()
	out_path = '/Users/Hugh/Documents/C1_test'
	ANIR,N = MyGetANIR(sta1,sta2,comp1,comp2,out_path)
	print ">>> Job Finished. It took %.2f seconds" %(time.time()-jt)
	#plt.plot(ANIR)
	#plt.show()

	pair = "%s:%s" % (sta1,sta2)
	out_path = out_path + '/'+comp1+'_'+comp2+'/'+sta1+'/'+'all_year/'
	filename = out_path + sta1 + '_' + sta2 + '.SAC'
	if N >0:
		export_sac(ANIR,filename,sta1,sta2,N)
		del ANIR,N
	#try:
        #	os.makedirs(os.path.split(filename)[0])
	#except:
        #	pass

	#mytrace = Trace(data=ANIR)
	#mytrace.stats['station'] = pair
	#mytrace.stats.sac = AttribDict()
	#mytrace.stats.sac.npts = len(ANIR)

	#st = Stream(traces=[mytrace,])
	#st.write(filename,format ='SAC')
	#del st


if __name__ == "__main__":
	jt= time.time()
	comp1 = 'BHZ'
	comp2 = 'BHZ'
	stations = model.load_stations('sources.txt')
	receivers = model.load_stations('receivers.txt')
	params = {}
   	params['sources_to_corr']       = ['%s.%s'%(sta.network,sta.station) for sta in stations]
	params['receivers_to_corr']     = ['%s.%s'%(sta.network,sta.station) for sta in receivers]
	pairs_to_compute = []
	for ia, a in enumerate(params['sources_to_corr']):
    		for b in params['receivers_to_corr']:
        		pairs_to_compute.append('%s:%s'%(a,b))
	print pairs_to_compute
	for i in range(len(pairs_to_compute)):
		sta1,sta2 = pairs_to_compute[i].split(':')
		run(sta1,sta2,comp1,comp2)
	
#	print ">>> all the jobs are finished. It took %.2f seconds in total" %(time.time()-jt)
