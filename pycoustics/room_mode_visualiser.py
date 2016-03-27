def calculate_modes(dimensions, maxfreq=2000):
    # to implement: arguments should be up to 3 dimensions in meters (check format and type)
    # extra: choose metres or feet
    # extra: temperature (influence on speed of sound)
    
    frequencies = [] # mode frequencies
    
    for arg in range(1,len(dimensions)):
        frequencies.append([]) # append empty list
        root = 170/float(dimensions[arg]) # f_0 = 340 m/s / (2*L)
        frequencies[arg-1].append(root)
        multiplier = 2
        while multiplier*root < maxfreq:
            frequencies[arg-1].append(root*multiplier)
            multiplier+=1 # to next harmonic
    
    return frequencies

def plot_room_modes(frequencies):
    # plot frequencies in red, blue and green, on logarithmic frequency axis
    # TODO extra: highlight region where near-coincidence occurs (show frequency)
    fig = plt.figure()
    colormap = ['b', 'r', 'g']
    for arg in range(1,len(sys.argv)):
        for line in range(0,len(frequencies[arg-1])):
            plt.plot([frequencies[arg-1][line], frequencies[arg-1][line]],
                     [.5*arg-1, .5*arg],
                     color=colormap[arg-1]
                     )
    
    plt.xlabel('Frequency (Hz)')
    plt.xscale('log')
    plt.title('Room modes')
    plt.xlim(10, maxfreq)
    
    # remove ticks
    frame = plt.gca()
    frame.axes.get_yaxis().set_visible(False)
    
    return fig
