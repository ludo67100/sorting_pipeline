#from compute_timefreq import compute_timefreq
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy import fftpack

#---------------------------------------!!!!!!!!!! FUNCTIONS !!!!!!!!!!-------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#---------------------------------------------------DO NOT--------------------------------------------------------
#---------------------------------------------------MODIFY--------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#--------------------------------------------GO BELOW FOR DATA INPUT----------------------------------------------
def generate_wavelet_fourier(len_wavelet, f_start, f_stop, delta_freq, 
            sampling_rate, f0, normalisation):
    """
    Compute the wavelet coefficients at all scales and makes its Fourier transform.
    When different signal scalograms are computed with the exact same coefficients, 
        this function can be executed only once and its result passed directly to compute_morlet_scalogram
        
    Output:
        wf : Fourier transform of the wavelet coefficients (after weighting), Fourier frequencies are the first 
    """
    # compute final map scales
    scales = f0/np.arange(f_start,f_stop,delta_freq)*sampling_rate
    # compute wavelet coeffs at all scales
    xi=np.arange(-len_wavelet/2.,len_wavelet/2.)
    xsd = xi[:,np.newaxis] / scales
    wavelet_coefs=np.exp(complex(1j)*2.*np.pi*f0*xsd)*np.exp(-np.power(xsd,2)/2.)

    weighting_function = lambda x: x**(-(1.0+normalisation))
    wavelet_coefs = wavelet_coefs*weighting_function(scales[np.newaxis,:])

    # Transform the wavelet into the Fourier domain
    #~ wf=fft(wavelet_coefs.conj(),axis=0) <- FALSE
    wf=fftpack.fft(wavelet_coefs,axis=0)
    wf=wf.conj() # at this point there was a mistake in the original script
    
    return wf


def convolve_scalogram(sig, wf):
    """
    Convolve with fft the signal (in time domain) with the wavelet
    already computed in freq domain.
    
    Parameters
    ----------
    sig: numpy.ndarray (1D, float)
        The signal
    wf: numpy.array (2D, complex)
        The wavelet coefficient in fourrier domain.
    """
    n = wf.shape[0]
    assert sig.shape[0]<=n, 'the sig.size is longer than wf.shape[0] {} {}'.format(sig.shape[0], wf.shape[0])
    sigf=fftpack.fft(sig,n)
    wt_tmp=fftpack.ifft(sigf[:,np.newaxis]*wf,axis=0)
    wt = fftpack.fftshift(wt_tmp,axes=[0])
    return wt

def compute_timefreq(sig, sampling_rate, f_start, f_stop, delta_freq=1., nb_freq=None,
                f0=2.5,  normalisation = 0.,  min_sampling_rate=None, wf=None,
                t_start=0., zero_pad=True, joblib_memory=None):
    """
    
    """
    #~ print 'compute_timefreq'
    sampling_rate = float(sampling_rate)
    
    if nb_freq is not None:
        delta_freq = (f_stop-f_start)/nb_freq
    
    if min_sampling_rate is None:
        min_sampling_rate =  min(4.* f_stop, sampling_rate)
        
    
    #decimate
    ratio = int(sampling_rate/min_sampling_rate)
    #~ print 'ratio', ratio
    if ratio>1:
        # sig2 = tools.decimate(sig, ratio)
        sig2 = scipy.signal.decimate(sig, ratio,n=4, zero_phase=True) #ORDER OF 4 FOR SMALL FREQ INTERVALLS !!!!!
    else:
        sig2 = sig
        ratio=1
    
    tfr_sampling_rate = sampling_rate/ratio
    #~ print 'tfr_sampling_rate', tfr_sampling_rate
    
    n_init = sig2.size
    if zero_pad:
        n = int(2 ** np.ceil(np.log(n_init)/np.log(2))) # extension to next power of  2
    else:
        n = n_init
    #~ print 'n_init', n_init, 'n', n
    if wf is None:
        if joblib_memory is None:
            func = generate_wavelet_fourier
        else:
            func = joblib_memory.cache(generate_wavelet_fourier)
        wf = func(n, f_start, f_stop, delta_freq, 
                            tfr_sampling_rate, f0, normalisation)
    
    assert wf.shape[0] == n
    
    wt = convolve_scalogram(sig2, wf)
    wt=wt[:n_init,:]
    
    freqs = np.arange(f_start,f_stop,delta_freq)
    times = np.arange(n_init)/tfr_sampling_rate + t_start
    return wt, times, freqs, tfr_sampling_rate


#-------------------------RIDGE EXTRACTION in 2D-------------------------------
#------------------------------------------------------------------------------
def ridge_map(ampl_map, threshold=70.):
    max_power = np.max(ampl_map) #Max power observed in frequency spectrum
    freq_power_threshold = float(threshold) #The threshold range for power detection of the ridge 
    cut_off_power = max_power/100.0*freq_power_threshold #Computes power above trheshold
    
    boolean_map = ampl_map >= cut_off_power #For plot 
    
    value_map = ampl_map
    
    for i,j in np.ndenumerate(ampl_map):
        if j <= cut_off_power:
            value_map[i] = 0.0 #For computation, all freq < trhesh = 0.0
            
    return boolean_map, value_map


#-------------------------SIGNAL FILTERING FUNCTION----------------------------
#------------------------------------------------------------------------------
def filter_signal(signal, order=8, sample_rate=20000,freq_low=400,freq_high=2000, axis=0):
    
    import scipy.signal
    
    Wn = [freq_low/(sample_rate/2),freq_high/(sample_rate/2)]
    
    sos_coeff = scipy.signal.iirfilter(order,Wn,btype='band',ftype='butter',output='sos')
    
    filtered_signal = scipy.signal.sosfiltfilt(sos_coeff, signal, axis=axis)
    
    return filtered_signal

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':  

#--------------------------------------------- Signal info : FILL HERE -------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
    #Location of the file
    file_loc ='H:/Federica/RAW_BINARY_FILES/2019-03-06T17-55-00McsRecording_1300um_P3.rbf'
    
    #Sampling rate (in Hertz) and start time (in seconds) and channel
    sampling_rate= 20000. 
    t_start = 0. 
    ch = 0
    
    #Freq window for Morlet Wavelet scalo (in Hz)
    f_start = 1.
    f_stop = 100.
    
    #Freq borders for filtered plot (in Hz)
    freq_low=600
    freq_high=2000
    
    #Treshold for ridge dection (% of the max power spectrum) 
    threshold = 50
    
#--------------------------------------------- THAT'S IT-------------- -------------------------------------------
#-----------------------------------------------------------------------------------------------------------------  
    #Define the signal, its length in seconds and a time vector
    sig = np.fromfile(file_loc,dtype='float64').reshape(-1,16)    
    duration = 1./sampling_rate * len(sig[:,ch])
    sig_times = np.arange(0, duration, 1./sampling_rate)
    
    
    #Fig 1 : the raw signal
    fig, ax = plt.subplots()
    ax.set_title('Here is the signal')
    ax.plot(sig_times, sig[:,ch], label='Raw signal')
    ax.plot(sig_times,filter_signal(sig[:,ch],freq_low=freq_low,freq_high=freq_high), label='Filtered signal')
    ax.legend(loc='best')
    
    #Fig 2 : Timefreq
    fig, ax = plt.subplots(1,2)
    ax[0].set_xlabel('Time (s)')
    ax[0].set_ylabel('Freq (Hz)')
    ax[0].set_title('Scalogram')
    
    complex_map, map_times, freqs, tfr_sampling_rate = compute_timefreq(sig[:,ch], sampling_rate, f_start, f_stop, delta_freq=1., nb_freq=None,
                    f0=2.5,  normalisation = 0., t_start=t_start)
    
    ampl_map = np.abs(complex_map) # the amplitude map (module)
    phase_map = np.angle(complex_map) # the phase
    
    delta_freq = freqs[1] - freqs[0]
    extent = (map_times[0], map_times[-1], freqs[0]-delta_freq/2., freqs[-1]+delta_freq/2.)
    
    scalo = ax[0].imshow(ampl_map.transpose(), interpolation='nearest', 
                        origin ='lower', aspect = 'auto', extent = extent, cmap='viridis')

    fig.colorbar(scalo)
    
    ridge_map_plot, ridge_map_comp = ridge_map(ampl_map, threshold=threshold)
    
    ax[1].set_title('Ridge freq ({} % of max power)'.format(threshold))
    ax[1].set_xlabel('Time (s)')
    ridge= ax[1].imshow(ridge_map_plot.transpose(), interpolation='nearest', 
                        origin ='lower', aspect = 'auto', extent = extent, cmap='hot')
    

    
    
    '''
    ind_max = np.argmax(ampl_map, axis=1)
    #print(ind_max)
    #~ phase_with_morlet = 
    freq_max = freqs[ind_max]
    ax.plot(map_times, freq_max, color='m', lw=1)
    '''

'''
fig, ax =  plt.subplots()
im = ax.imshow(phase_map.transpose(), interpolation='nearest', 
                    origin ='lower', aspect = 'auto', extent = extent, cmap='viridis')
im.set_clim(-np.pi, np.pi)
fig.colorbar(im)
ax.plot(map_times, freq_max, color='m', lw=3)



estimated_phase_morlet = phase_map[np.arange(map_times.size), ind_max]
estimated_ampl_morlet = ampl_map[np.arange(map_times.size), ind_max]

# hilbert way
analytic_signal = scipy.signal.hilbert(sig)
amplitude_envelope = np.abs(analytic_signal)
instantaneous_phase = np.unwrap(np.angle(analytic_signal))
instantaneous_frequency = (np.diff(instantaneous_phase) / (2.0*np.pi) * sampling_rate)



fig, axs = plt.subplots(nrows=2)
axs[0].plot(sig_times, sig_ampl, color='g')
axs[0].plot(map_times, estimated_ampl_morlet, color='r')
axs[0].plot(sig_times, amplitude_envelope, color='c')

axs[1].plot(sig_times, sig_phase%(2*np.pi), color='g')
axs[1].plot(map_times,estimated_phase_morlet, color='r')
axs[1].plot(sig_times,instantaneous_phase%(2*np.pi), color='c')


##### spike phase
spike_indexes_400 = (spike_times*tfr_sampling_rate).astype('int64')

spikes_phase_morlet = estimated_phase_morlet[spike_indexes_400]
spikes_phase_hilbert = instantaneous_phase[spike_indexes]

spikes_phase_morlet = (spikes_phase_morlet+2*np.pi) % (2*np.pi)
spikes_phase_hilbert = (spikes_phase_hilbert+2*np.pi) % (np.pi*2)

bins=np.arange(-np.pi, 3*np.pi, np.pi/10)
hist_phase_morlet, bins = np.histogram(spikes_phase_morlet, bins=bins)
hist_phase_hilbert, bins = np.histogram(spikes_phase_hilbert, bins=bins)


fig, ax = plt.subplots()
ax.plot(bins[:-1], hist_phase_morlet, color='r')
ax.plot(bins[:-1], hist_phase_hilbert, color='c')
'''





