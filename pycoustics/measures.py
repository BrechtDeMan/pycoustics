def EDC(impulse_response, fs):
    fs = float(fs) # make sure this is treated as a floating point number
    impulse_response = np.array(impulse_response) # convert to np array if it isn't already
    test_monophonic(impulse_response)
    
    # Schroeder's energy decay curve
    cumul = 10.0*np.log10(np.sum(np.square(impulse_response)))
    decay_curve = 10.0*np.log10(np.flipud(np.cumsum(np.flipud(np.square(impulse_response))))) - cumul
    return decay_curve
    
def RT(decay_curve, fs, end_level, start_level=-5):
    t_start = np.argmax(decay_curve<start_level)/fs # time at which EDC drops below -5dB or start_level
    t_end   = np.argmax(decay_curve<end_level)/fs # time at which EDC drops below end_level
    if t_5==0 or t_25==0: # argmax==0 when condition never satisfied
        RT_d=RT_r=float('nan')
    else: 
        # decay method
        RT_d = (t_end-t_start)*60/(end_level-start_level)
        # regression method
        s_start = np.ceil(t_start*fs).astype(int) # convert start time to integer sample number
        s_end = np.ceil(t_end*fs).astype(int) # convert end time to integer sample number
        p = np.polyfit(np.arange(s_start,s_end), decay_curve[np.arange(s_start,s_end)],1)
        RT_r = -60/p[0]/fs
    return (RT_d, RT_r, p)

def RT20(decay_curve, fs):
    return RT(decay_curve, fs, -25)
    
def RT30(decay_curve, fs):
    return RT(decay_curve, fs, -35)
    
def EDT(decay_curve, fs):
    return RT(decay_curve, fs, -10, 0)

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
