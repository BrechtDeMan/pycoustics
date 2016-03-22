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

def K_filter(signal, fs, debug=False):
    # apply K filtering as specified in EBU R-128 / ITU BS.1770-4

    # pre-filter 1
    f0 = 1681.9744509555319
    G  = 3.99984385397
    Q  = 0.7071752369554193
    K  = np.tan(np.pi * f0 / fs) # TODO: precompute
    Vh = np.power(10.0, G / 20.0)
    Vb = np.power(Vh, 0.499666774155)
    a0_ = 1.0 + K / Q + K * K
    b0 = (Vh + Vb * K / Q + K * K) / a0_
    b1 = 2.0 * (K * K -  Vh) / a0_
    b2 = (Vh - Vb * K / Q + K * K) / a0_
    a0 = 1.0
    a1 = 2.0 * (K * K - 1.0) / a0_
    a2 = (1.0 - K / Q + K * K) / a0_
    signal_1 = lfilter([b0,b1,b2],[a0,a1,a2],signal)
    
    if debug:
        plt.figure(figsize=(9,9))
        ax1 = fig.add_subplot(111)
        w, h1 = freqz([b0,b1,b2], [a0,a1,a2], worN=8000)#np.logspace(-4, 3, 2000))
        plt.semilogx((fs * 0.5 / np.pi) * w, 20*np.log10(abs(h1)))
        plt.title('Pre-filter 1')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Gain [dB]')
        plt.xlim([20, 20000])
        plt.ylim([-10,10])
        plt.grid(True, which='both')
        ax = plt.axes()
        ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        plt.show()
    
    # pre-filter 2
    f0 = 38.13547087613982
    Q  =  0.5003270373253953
    K  = np.tan(np.pi * f0 / fs)
    a0 = 1.0
    a1 = 2.0 * (K * K - 1.0) / (1.0 + K / Q + K * K)
    a2 = (1.0 - K / Q + K * K) / (1.0 + K / Q + K * K)
    b0 = 1.0
    b1 = -2.0
    b2 = 1.0
    signal_2 = lfilter([b0,b1,b2],[a0,a1,a2],signal_1)
    
    if debug:
        plt.figure(figsize=(9,9))
        ax1 = fig.add_subplot(111)
        w, h2 = freqz([b0,b1,b2], [a0,a1,a2], worN=8000)#np.logspace(-4, 3, 2000))
        plt.semilogx((fs * 0.5 / np.pi) * w, 20*np.log10(abs(h2)))
        plt.title('Pre-filter 2')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Gain [dB]')
        plt.xlim([10, 20000])
        plt.ylim([-30,5])
        plt.grid(True, which='both')
        ax = plt.axes()
        ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
        plt.show()
    
    return signal_2 # return signal passed through 2 pre-filters
