def C50(impulse_response, fs, delay=0):
    return clarity(.050, impulse_response, fs, delay)

def C80(impulse_response, fs, delay=0):
    return clarity(.080, impulse_response, fs, delay)

def clarity(early_time, impulse_response, fs, delay=0):
    D = definition(early_time, impulse_response, fs, delay)
    return 10.0*np.log10(D/(1-D))

def D50(impulse_response, fs, delay=0):
    return definition(.050, impulse_response, fs, delay)

def definition(early_time, impulse_response, fs, delay=0):
    test_monophonic(impulse_response) # how to call this function from another file? 
    # disregard everything before delay
    impulse_response = impulse_response[np.round(delay*fs).astype(int):] 
    early_energy = np.sum(np.square(impulse_response[:np.round(early_time*fs).astype(int)]))
    total_energy = np.sum(np.square(impulse_response))
    return early_energy/total_energy

def TS(impulse_response, fs, delay=0): # centre time
    test_monophonic(impulse_response)
    # disregard everything before delay
    impulse_response = impulse_response[np.round(delay*fs).astype(int):] 
    numerator = np.sum(np.arange(len(impulse_response))/fs * np.square(impulse_response)) # ∑ t * p^2(t)
    total_energy = np.sum(np.square(impulse_response)) # ∑ p^2(t)
    return 1000.0*numerator/total_energy # expressed in ms
    
def calculate_loudness(signal,fs):
    # filter
    if len(signal.shape)==1: # if shape (N,), then make (N,1)
        signal_filtered = copy.copy(signal.reshape((signal.shape[0],1)))
    else:
        signal_filtered = copy.copy(signal)
        
    for i in range(signal_filtered.shape[1]):
        signal_filtered[:,i] = K_filter(signal_filtered[:,i], fs)

    # mean square
    G = [1.0, 1.0, 1.0, 1.41, 1.41]
    T_g = 0.400 # 400 ms gating block
    Gamma_a = -70.0 # absolute threshold: -70 LKFS
    overlap = .75 # relative overlap (0.0-1.0)
    step = 1 - overlap

    T = signal_filtered.shape[0]/fs # length of measurement interval in seconds
    j_range = np.arange(0,(T-T_g)/(T_g*step))
    z = np.ndarray(shape=(signal_filtered.shape[1],len(j_range))) # ?
    # write in explicit for-loops for readability and translatability
    for i in range(signal_filtered.shape[1]): # for each channel i
        for j in j_range: # for each window j
            lbound = np.round(fs*T_g*j*step).astype(int)
            hbound = np.round(fs*T_g*(j*step+1)).astype(int)
            z[i,j] = (1/(T_g*fs))*np.sum(np.square(signal_filtered[lbound:hbound, i]))

    G_current = np.array(G[:signal_filtered.shape[1]]) # discard weighting coefficients G_i unused channels
    n_channels = G_current.shape[0]
    l = [-.691 + 10.0*np.log10(np.sum([G_current[i]*z[i,j.astype(int)] for i in range(n_channels)])) \
             for j in j_range]
    #print 'l: '+str(l)

    # throw out anything below absolute threshold:
    indices_gated = [idx for idx,el in enumerate(l) if el > Gamma_a] 
    z_avg = [np.mean([z[i,j] for j in indices_gated]) for i in range(n_channels)]
    Gamma_r = -.691 + 10.0*np.log10(np.sum([G_current[i]*z_avg[i] for i in range(n_channels)])) - 10.0
    # throw out anything below relative threshold:
    indices_gated = [idx for idx,el in enumerate(l) if el > Gamma_r] 
    z_avg = [np.mean([z[i,j] for j in indices_gated]) for i in range(n_channels)]
    L_KG = -.691 + 10.0*np.log10(np.sum([G_current[i]*z_avg[i] for i in range(n_channels)]))

    return L_KG

def RT(decay_curve, fs, end_level, start_level=-5):
    t_start = np.argmax(decay_curve<start_level)/fs # time at which EDC drops below -5dB or start_level
    t_end   = np.argmax(decay_curve<end_level)/fs # time at which EDC drops below end_level
    if t_5==0 or t_25==0: # argmax==0 when condition never satisfied
        RT_d=RT_r=float('nan')
        p2=[float('nan'), float('nan')]
    else: 
        # decay method
        RT_d = (t_end-t_start)*60/(end_level-start_level)
        # regression method
        s_start = np.ceil(t_start*fs).astype(int) # convert start time to integer sample number
        s_end = np.ceil(t_end*fs).astype(int) # convert end time to integer sample number
        p = np.polyfit(np.arange(s_start,s_end), decay_curve[np.arange(s_start,s_end)],1)
        RT_r = -60/p[0]/fs
    return (RT_d, RT_r)
