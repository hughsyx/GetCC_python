from GetANIR import *
import time
from obspy.core import Trace, Stream, AttribDict
from pyrocko import model
import pdb
from multiprocessing import Process

def run(pair,comp1,comp2):
	jt = time.time()
	sta1,sta2 = pair.split(':')
	out_path = '../../C1_SC'
	MyGetANIR(sta1,sta2,comp1,comp2,out_path)
	print ">>> Job Finished. It took %.2f seconds for %s" %(time.time()-jt,pair)
	#plt.plot(ANIR)
	#plt.show()

#	pair = "%s:%s" % (sta1,sta2)
#	out_path = '../../C1_SC/'+comp1+'_'+comp2+'/'+sta1+'/'
#	filename = out_path + sta1 + '_' + sta2 + '.SAC'
#	try:
#        	os.makedirs(os.path.split(filename)[0])
#	except:
#        	pass

#	mytrace = Trace(data=ANIR)
#	mytrace.stats['station'] = pair
#	mytrace.stats.sac = AttribDict()
#	mytrace.stats.sac.npts = len(ANIR)

#	st = Stream(traces=[mytrace,])
#	st.write(filename,format ='SAC')
#	del st


if __name__ == "__main__":
	jt_all = time.time()
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
	nthreads = 12
	clients = []
	for pair in pairs_to_compute:
		client = Process(target=run,args=(pair,comp1,comp2))
		clients.append(client)
		client.start()
		
		while len(clients) >= nthreads:
			for client in clients:
				client.join(0.01)
				if not client.is_alive():
					client.join(0.01)
					clients.remove(client)
	while len(clients) != 0:
		for client in clients:
			client.join(0.01)
			if not client.is_alive():
				client.join()
				clients.remove(client)

	print ">>> All the jobs have been finished. It took %.2f seconds in total" %(time.time()-jt_all)
