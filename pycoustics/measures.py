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
