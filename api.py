import os
import copy
import datetime
import itertools
import math
import pandas as pd

import numpy as np
import scipy.fftpack
from obspy.core import Stream, Trace, read, AttribDict
from obspy.signal.invsim import cosTaper

from obspy.core.util import gps2DistAzimuth
from pyrocko import model


def datetime_range(start, end, delta):
    current = start
    if not isinstance(delta, datetime.timedelta):
        delta = datetime.timedelta(**delta)
    while current < end:
        yield current
        current += delta

def moving_avg(x,n):
    n = int(n)
    cx = x.cumsum()
    nn = len(x)
    y = np.zeros(nn, dtype=cx.dtype)
    y[n/2:n/2+(nn-n)] = (cx[n:]-cx[:-n])/n
    y[:n/2] = y[n/2]
    y[n/2+(nn-n):] = y[n/2+(nn-n)-1]
    return y


def get_station(sta):
    stations = model.load_stations('stations.txt')
    for s in stations:
        if s.station == sta:
            return s
            break



def get_interstation_distance(station1, station2, coordinates="DEG"):
    """Returns the distance in km between `station1` and `station2`.

    .. warning:: Currently the stations coordinates system have to be the same!

    :type station1: :class:`~msnoise.msnoise_table_def.Station`
    :param station1: A Station object
    :type station2: :class:`~msnoise.msnoise_table_def.Station`
    :param station2: A Station object
    :type coordinates: str
    :param coordinates: The coordinates system. "DEG" is WGS84 latitude/
        longitude in degrees. "UTM" is expressed in meters.



    :rtype: float
    :returns: The interstation distance in km
    """

    if coordinates == "DEG":
        dist, azim, bazim = gps2DistAzimuth(station1.Y, station1.X,
                                            station2.Y, station2.X)
        return dist / 1.e3
    else:
        dist = np.hypot(float(station1.X - station2.X),
                        float(station1.Y - station2.Y)) / 1.e3
        return dist


# CORRELATIONS

#***
def export_allcorr(params, ccfid, data, name='corr'):
    output_folder = params['output_folder']
    station1, station2, filterid, components, date = ccfid.split('_')

    path = os.path.join(output_folder,'%s'%name, "%02i" % int(filterid),
                        station1, station2, components)
    if not os.path.isdir(path):
        os.makedirs(path)

    df = pd.DataFrame().from_dict(data).T
    print df
    df.columns = get_t_axis(params)
    df.to_hdf(os.path.join(path, date+'.h5'), 'data')
    del df
    return


def add_corr(params, station1, station2, filterid, date, time, duration, components, CF, sampling_rate,name='corr', day=False, ncorr=0):
    """
    Adds a CCF to the data archive on disk.
    
    :type params: dict 
    :param params: This dictionary contains all parameters for correlaion. see params.py to initilalize it.
    :type station1: str
    :param station1: The name of station 1 (formatted NET.STA)
    :type station2: str
    :param station2: The name of station 2 (formatted NET.STA)
    :type filterid: int
    :param filterid: The ID (ref) of the filter
    :type date: datetime.date or str
    :param date: The date of the CCF
    :type time: datetime.time or str
    :param time: The time of the CCF
    :type duration: float
    :param duration: The total duration of the exported CCF
    :type components: str
    :param components: The name of the components used (ZZ, ZR, ...)
    :type sampling_rate: float
    :param sampling_rate: The sampling rate of the exported CCF
    :type day: bool
    :param day: Whether this function is called to export a daily stack (True)
        or each CCF (when keep_all parameter is set to True in the
        configuration). Defaults to True.
    :type ncorr: int
    :param ncorr: Number of CCF that have been stacked for this CCF.
    """

    output_folder = params['output_folder']
    export_format = params['export_format']
    sac, mseed = False, False
    if export_format in ["BOTH","both"]:
        mseed = True
        sac = True
    elif export_format in ["SAC","sac"]:
        sac = True
    elif export_format in ["MSEED","mseed"]:
        mseed = True
    if params['crosscorr']:pass
    
    if day:
        path = os.path.join("STACKS",'%s'%name , "%02i" % filterid, "001_DAYS", components,
                            "%s_%s" % (station1, station2), str(date))
        pair = "%s:%s" % (station1, station2)
        if mseed:
            export_mseed(params, path, pair, components, filterid, CF/ncorr,
                         ncorr)
        if sac:
            export_sac(params, path, pair, components, filterid, CF/ncorr,
                       ncorr)

    else:
        file = '%s.cc' % time
        path = os.path.join(output_folder, "%02i" % filterid, station1,
                            station2, components, date)
        if not os.path.isdir(path):
            os.makedirs(path)

        t = Trace()
        t.data = CF
        t.stats.sampling_rate = sampling_rate
        t.stats.starttime = -float(params['maxlag'])
        t.stats.components = components
        # if ncorr != 0:
            # t.stats.location = "%02i"%ncorr
        st = Stream(traces=[t, ])
        st.write(os.path.join(path, file), format='mseed')
        del t, st


def export_sac(params, filename, pair, components, filterid, corr, ncorr=0,
               sac_format=None, maxlag=None, cc_sampling_rate=None):
    if sac_format is None:
        sac_format = params["sac_format"]
    if maxlag is None:
        maxlag = float(params["maxlag"])
    if cc_sampling_rate is None:
        cc_sampling_rate = float(params["cc_sampling_rate"])
    try:
        os.makedirs(os.path.split(filename)[0])
    except:
        pass
    filename += ".sac"
    mytrace = Trace(data=corr)
    mytrace.stats['station'] = pair
    mytrace.stats['sampling_rate'] = cc_sampling_rate
    if maxlag:
        mytrace.stats.starttime = -maxlag
    mytrace.stats.sac = AttribDict()
    mytrace.stats.sac.depmin = np.min(corr)
    mytrace.stats.sac.depmax = np.max(corr)
    mytrace.stats.sac.depmen = np.mean(corr)
    mytrace.stats.sac.scale = 1
    mytrace.stats.sac.npts = len(corr)

    st = Stream(traces=[mytrace, ])
    st.write(filename, format='SAC')
    del st
    return


def export_mseed(params, filename, pair, components, filterid, corr, ncorr=0,
                 maxlag=None, cc_sampling_rate=None):
    try:
        os.makedirs(os.path.split(filename)[0])
    except:
        pass
    filename += ".mseed"
    
    if maxlag is None:
        maxlag = float(params["maxlag"])
    if cc_sampling_rate is None:
        cc_sampling_rate = float(params["cc_sampling_rate"])

    mytrace = Trace(data=corr)
    mytrace.stats['station'] = pair[:11]
    mytrace.stats['sampling_rate'] = cc_sampling_rate
    mytrace.stats['start_time'] = -maxlag
    mytrace.stats['location'] = "%02i" % ncorr

    st = Stream(traces=[mytrace, ])
    st.write(filename, format='MSEED')
    del st
    return


# Some helper functions

#*******
def get_maxlag_samples(params):
    """
    Returns the length of the CC functions. Gets the maxlag and sampling rate
    from the database.


    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: int
    :returns: the length of the CCF
    """

    maxlag = float(params['maxlag'])
    cc_sampling_rate = float(params['cc_sampling_rate'])

    if maxlag == params['corr_duration']:
        return int(2*maxlag*cc_sampling_rate)-1
    else: 
        return int(2*maxlag*cc_sampling_rate)+1

#****
def get_t_axis(params):
    """
    Returns the time axis (in seconds) of the CC functions.
    Gets the maxlag from the database and uses `get_maxlag_samples` function.

    :rtype: :class:`numpy.array`
    :returns: the time axis
    """

    maxlag = float(params['maxlag'])
    samples = get_maxlag_samples(params)
    return np.linspace(-maxlag, maxlag, samples)



# MISC


def azimuth(x0, y0, x1, y1):
    """
    Returns the azimuth between two coordinate sets.

    :type x0: float
    :param x0: X coordinate of station 1
    :type y0: float
    :param y0: Y coordinate of station 1
    :type x1: float
    :param x1: X coordinate of station 2
    :type y1: float
    :param y1: Y coordinate of station 2

    :rtype: float
    :returns: The azimuth in degrees
    """
    
    dist, azim, bazim = gps2DistAzimuth(y0, x0, y1, x1)
    return azim


def nextpow2(x):
    """
    Returns the next power of 2 of `x`.

    :type x: int
    :param x: any value

    :rtype: int
    :returns: the next power of 2 of `x`
    """

    return np.ceil(np.log2(np.abs(x)))


def check_and_phase_shift(trace):
    # print trace
    taper_length = 20.0
    if trace.stats.npts < 4 * taper_length*trace.stats.sampling_rate:
        trace.data = np.zeros(trace.stats.npts)
        return trace

    dt = np.mod(trace.stats.starttime.datetime.microsecond*1.0e-6,
                trace.stats.delta)
    if (trace.stats.delta - dt) <= np.finfo(float).eps:
        dt = 0
    if dt != 0:
        if dt <= (trace.stats.delta / 2.):
            dt = -dt
#            direction = "left"
        else:
            dt = (trace.stats.delta - dt)
#            direction = "right"
        trace.detrend(type="demean")
        trace.detrend(type="simple")
        taper_1s = taper_length * float(trace.stats.sampling_rate) / trace.stats.npts
        cp = cosTaper(trace.stats.npts, taper_1s)
        trace.data *= cp

        n = int(2**nextpow2(len(trace.data)))
        FFTdata = scipy.fftpack.fft(trace.data, n=n)
        fftfreq = scipy.fftpack.fftfreq(n, d=trace.stats.delta)
        FFTdata = FFTdata * np.exp(1j * 2. * np.pi * fftfreq * dt)
        trace.data = np.real(scipy.fftpack.ifft(FFTdata, n=n)[:len(trace.data)])
        trace.stats.starttime += dt
        return trace
    else:
        return trace


def getGaps(stream, min_gap=None, max_gap=None):
    # Create shallow copy of the traces to be able to sort them later on.
    copied_traces = copy.copy(stream.traces)
    stream.sort()
    gap_list = []
    for _i in xrange(len(stream.traces) - 1):
        # skip traces with different network, station, location or channel
        if stream.traces[_i].id != stream.traces[_i + 1].id:
            continue
        # different sampling rates should always result in a gap or overlap
        if stream.traces[_i].stats.delta == stream.traces[_i + 1].stats.delta:
            flag = True
        else:
            flag = False
        stats = stream.traces[_i].stats
        stime = stats['endtime']
        etime = stream.traces[_i + 1].stats['starttime']
        delta = etime.timestamp - stime.timestamp
        # Check that any overlap is not larger than the trace coverage
        if delta < 0:
            temp = stream.traces[_i + 1].stats['endtime'].timestamp - \
                etime.timestamp
            if (delta * -1) > temp:
                delta = -1 * temp
        # Check gap/overlap criteria
        if min_gap and delta < min_gap:
            continue
        if max_gap and delta > max_gap:
            continue
        # Number of missing samples
        nsamples = int(round(math.fabs(delta) * stats['sampling_rate']))
        # skip if is equal to delta (1 / sampling rate)
        if flag and nsamples == 1:
            continue
        elif delta > 0:
            nsamples -= 1
        else:
            nsamples += 1
        gap_list.append([_i, _i+1,
                        stats['network'], stats['station'],
                        stats['location'], stats['channel'],
                        stime, etime, delta, nsamples])
    # Set the original traces to not alter the stream object.
    stream.traces = copied_traces
    return gap_list


if __name__ == "__main__":
    pass

#EOF
