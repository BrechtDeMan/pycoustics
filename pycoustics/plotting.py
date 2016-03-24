def plot_impulse_response(impulse_response, fs, win_length=0): 
    irsq = np.square(impulse_response) # squared impulse response
    if win_length > 0:
        window = np.ones(np.round(win_length*fs).astype(int))
        irsq = np.convolve(irsq, window) # smoothing
    fig = plt.figure()
    plt.plot(np.arange(len(irsq))/fs, irsq/max(abs(irsq))) # plot impulse response
    plt.yscale('log')
    plt.xlabel('Time [s]')
    plt.ylabel('Amplitude [dB]')
    return fig
