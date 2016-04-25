from scipy.signal import butter, lfilter 

def butter_bandpass(lower_freq, higher_freq, fs, order=2):
    # based on http://scipy.github.io/old-wiki/pages/Cookbook/ButterworthBandpass
    nyq = 0.5 * fs # Nyquist frequency = half of sample rate
    b, a = butter(order, [lower_freq / nyq, higher_freq / nyq], btype='band')
    return b, a

def butter_cascade(x, lower_freq, higher_freq, fs, order=2):
    # filters signal with cascade of 2nd and 1st order Butterworth filters
    if order % 2 == 1: # if odd
        b1, a1 = butter_bandpass(lower_freq, higher_freq, fs, order=1)
        x = lfilter(b, a, x)
    for idx in range(order/2): # integer division
        b, a = butter_bandpass(lower_freq, higher_freq, fs, order=2)
        x = lfilter(b, a, x)
    return x

def octave_filter(x, fc, fs, order=2):
    sqrt2 = np.sqrt(2)
    lower_freq = fc/sqrt2
    higher_freq = fc*sqrt2
    y = butter_cascade(x, lower_freq, higher_freq, fs, order)
    return y

def third_octave_filter(x, fc, fs, order=2):
    sixth_root_of_2 = np.power(2.0,1.0/6)
    lower_freq = fc/sixth_root_of_2
    higher_freq = fc*sixth_root_of_2
    y = butter_cascade(x, lower_freq, higher_freq, fs, order)
    return y

def spectral_smoothing(signal, alpha): 
    # takes time domain signal and applies complex smoothing in spectral domain
    # simple: 
    signal_spectrum = np.fft.rfft(signal)
    nfft = len(signal_spectrum)
    if alpha>0.0:
        smoothing_window = np.hanning(alpha*nfft)
        signal_spectrum = np.convolve(signal_spectrum, smoothing_window, mode='same')
    # Maamar ICASSP (iterative): 
    #TODO ...
    return np.fft.irfft(signal_spectrum)
