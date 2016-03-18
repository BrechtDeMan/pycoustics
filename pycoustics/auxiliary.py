def seconds2timestr(time_in_seconds):
# Turn number of seconds (int) to '[minutes] min [seconds] s' (string)
    if time_in_seconds is not None and not isNaN(time_in_seconds): 
        time_in_minutes = int(time_in_seconds/60)
        remaining_seconds = int(time_in_seconds%60)
        return str(time_in_minutes) + " min " + str(remaining_seconds) + " s"
    else:
        return 'N/A'

def test_monophonic(signal):
    signal = np.array(signal) # in case not np.array
    assert signal.ndim == 1, 'Not a one dimensional vector.'
