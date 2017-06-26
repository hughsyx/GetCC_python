
from pyrocko import model
import logging

#Data Structure Definitions:
    #data_structure = {}
    #data_structure['SDS'] = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
    #data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
    #data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
    #data_structure['PDF'] = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"


def get_params():
    stations = model.load_stations('sources.txt')
    receivers = model.load_stations('receivers.txt')
    # Get Configuration
    params = {}
#    params['path_to_data']          = '/Volumes/VOLC_DATA/LASTARRIA/raw'#'/Users/Zack/Desktop/corrpy/RAW'
    params['path_to_data']          = '/Users/Hugh/Documents/data/SoCal'
    params['output_folder']         = '/Users/Hugh/Documents/CORR'
    params['architecture']          = 'IDDS'
    params['cc_sampling_rate']      = 5 #float(get_config(db, "cc_sampling_rate"))
    params['analysis_duration']     = 86400 # how much time to process as bunch. must stay like that not implemented to change now
    params['overlap']               = 0.5 
    params['maxlag']                = 1000 # even smaller window if needed
    params['corr_duration']         = 1800 # slicing the 86400 in small windows
    params['npts']                  = params['corr_duration'] * params['cc_sampling_rate']
    params['temp_norm']             = 0 # 0: removing eq with param 'clipping'; -1: 1-bit normalization; 1: Windsorinzing with param 'clipping'   
    params['clipping']              = 8 # clipping eq or windsorizing at 3 * std(trace)
    params['resampling_method']     = "Resample" #"Decimate"
    params['decimation_factor']     = int(5)
    params['preprocess_lowpass']    = 2.00#float(get_config(db, "preprocess_lowpass"))
    params['preprocess_highpass']   = 0.05#float(get_config(db, "preprocess_highpass"))
    params['keep_all']              = False#get_config(db, 'keep_all', isbool=True)
    params['keep_days']             = True#get_config(db, 'keep_days', isbool=True)
    params['components_to_compute'] = ['ZZ']#get_components_to_compute(db)
    params['sources_to_corr']       = ['%s.%s'%(sta.network,sta.station) for sta in stations]
    params['receivers_to_corr']     = ['%s.%s'%(sta.network,sta.station) for sta in receivers]
    params['starttime']             = '2014-01-01'
    params['endtime']               = '2014-12-31'
    params['export_format']         = 'sac'
    params['sac_format']            = 'doublets' #Format for SAC stacks ? [doublets]/clarke
    params['crosscorr']             = True
    params['deconvolution']         = True
    params['cross-coherence']       = False#True
    params['nthreads']              = 3

    filter1 = {} 
    filter1['ref']          = 1
    filter1['low']          = 0.05               # The lower frequency bound of the Whiten function (in Hz) 
    filter1['high']         = 2.00                  # The upper frequency bound of the Whiten function (in Hz)
    filter1['rms_threshold'] = 3
#    filter2 = {} 
#    filter2['ref']      = 2
#    filter2['low']      = 0.01               # The lower frequency bound of the Whiten function (in Hz) 
#    filter2['high']     = 8                  # The upper frequency bound of the Whiten function (in Hz)

    filters = {}
    filters['1']=filter1
#    filters['2']=filter2
    
    return params, filters


        



if __name__=='__main__':
    params, filters = get_params()
    print params
    print filters

#EDF
