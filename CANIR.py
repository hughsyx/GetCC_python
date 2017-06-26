import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
from obspy.signal.invsim import waterlevel
import pdb
import math
#from api import nextpow2

def myCorr(vb,va,npt):
	""" This function is used to get the cross correlation with deconvolution
	vb is the receiver and va is the virtual source
	npt is the number used for running window average
	"""
	eps = 1e-9
	lena = len(va)
	lenb = len(vb)
	Len = lena + lenb -1
	fva = scipy.fftpack.fft(va,Len)
	fvb = scipy.fftpack.fft(vb,Len)
	
	#fva2 = RWA2(fva,npt)
	#fvb2 = RWA2(fvb,npt)
	#p1 = np.ones(len(fva))*eps
	#p2 = np.ones(len(fvb))*eps
	
	#for i in range(Len):
	#	p1[i] = sum(fva2[i:i+npt])/npt;
	#	p2[i] = sum(fvb2[i:i+npt])/npt;
	
	#ANIRWF = fvb * np.conj(fva)/(p1**2)
	#ANIRW = np.real(scipy.fftpack.ifftshift(scipy.fftpack.ifft(ANIRWF)));
	#del ANIRWF, fva , fvb
	#return ANIRW
	p1 = moving_avg(abs(fva),npt)
	pdb.set_trace()
	ANIRWF = fvb * np.conj(fva)/(p1**2)
	ANIRW = np.real(scipy.fftpack.fftshift(scipy.fftpack.ifft(ANIRWF)));
	del ANIRWF, fva , fvb
	return ANIRW



def moving_avg(x,n):
        n = int(n)
        cx = x.cumsum()
        nn = len(x)
        y = np.zeros(nn, dtype=cx.dtype)
        y[n/2:n/2+(nn-n)] = (cx[n:]-cx[:-n])/n
        y[:n/2] = y[n/2]
        y[n/2+(nn-n):] = y[n/2+(nn-n)-1]
        return y


def myCorr2(vb,va,npt):
        """ This function is used to get the cross correlation with deconvolution
        vb is the receiver and va is the virtual source
        npt is the number used for running window average
    use power(2,N) as the length of fft
        """
        eps = 1e-9
        lena = len(va)
        lenb = len(vb)
        Len = lena+lenb-1
        Len = int(math.pow(2,nextpow2(Len)))
        fva = scipy.fftpack.fft(va,Len)
        fvb = scipy.fftpack.fft(vb,Len)

        p1 = moving_avg(abs(fva),npt)
        ANIRWF = fvb * np.conj(fva)/(p1**2)
        ANIRW = np.real(scipy.fftpack.fftshift(scipy.fftpack.ifft(ANIRWF)));
        ANIRW = ANIRW[Len/2-lena:Len/2+lenb-1]
        #ANIRW = ANIR[Len/2-lena+1:len/2+lenb]
	del ANIRWF, fva , fvb
        return ANIRW


def nextpow2(x):
    """
    Returns the next power of 2 of `x`.

    :type x: int
    :param x: any value

    :rtype: int
    :returns: the next power of 2 of `x`
    """

    return np.ceil(np.log2(np.abs(x)))



#def RWA2(v,npt):
#	""" this function is used to get the running window average
#	"""
#	Len = len(v)
#	v2 = np.zeros(Len + 2*npt) + np.zeros(Len + 2*npt)*0j
#	vtmp1 = v[0:npt]
#	vtmp2 = v[Len-npt:]
#	v2[0:npt] = vtmp1[::-1]
#	v2[npt+Len:] = vtmp2[::-1]
#	v2[npt:npt+Len]=v
#	v2 = abs(v2)
	
#	del v, vtmp1, vtmp2, Len
#	return v2

if __name__ == "__main__":
	
	data1 = np.random.random((1000,))
	data2 = np.random.random((1000,))
	#data1 = np.array([1,2,3,2,1,2,2,2,2,2,2,2,2,3,1,4,2,4,1,5,2,5,5,2,4,2,4,6,7,4,2,3,5])
	#data2 = np.array([2,3,4,1,2,4,1,4,1,3,5,2,2,5,1,3,5,1,3,5,3,5,4,2,3,6,3,3,6,7,3,4,6])
	corr = myCorr(data2,data1,10)
	print np.mean(corr)
	plt.figure()
	plt.plot(corr)
	plt.show()

