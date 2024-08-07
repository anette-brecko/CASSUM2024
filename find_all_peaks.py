import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
import pandas as pd
import time
start_time = time.time()

start_dir = '/Users/elsabrecko/Desktop/CASSUM/'
end_file = '_K.dat'

source_info = pd.read_csv(start_dir + 'source_info_converted.txt')
col = source_info.iloc[:, 0]
group_numbers = source_info.iloc[:, 5]

all_file_names = []
final_peak_per_GHz = []

for j in range(0, len(col), 4):
    #all_kelvin has all temperature values for the peaks
    all_kelvin = []
    
    #all_heights has the heights of the peaks in the order of which they're found
    all_heights = []
    
    #all_peaks has the frequencies of all the peaks
    all_peaks = []
    
    #There's one RMS value per SPW, all_rms keeps track of which RMS belongs to which source/SPW
    all_rms = []
    
    min_freq = float('inf')
    max_freq = float('-inf')
    
    #number_of_freq keeps track of how many frequencies are covered (minus overlap) in order to calcualte average peaks/GHz
    number_of_freq = 0
    
    #filename needs to include source number and MM number, which can be read in from source_info_converted.txt
    #Does not include SPW (inner for loop does this), so every 4th source name must be read from source_info_converted.txt
    filename = col[j][:-2]
    all_file_names.append(filename)

    #This inner for loop goes through all of the spectral windows, and allows for the min/max freq to be reset after a particular source
    for i in range(4):
        name = filename + '_' + str(i)
        RMS = source_info.loc[col == name, 'rms (Kelvin)'].values[0]
        dir = start_dir + 'SPARKS_spectra/GROUP' + str(group_numbers[j + i]) + '_K/'
        
        current_spw = pd.read_csv(dir + name + end_file)
        frequencies = list(current_spw.iloc[:, 0])
        kelvin_values = list(current_spw.iloc[:, 1])

        #Get index of frequencies that have not yet been covered by a different SPW
        #Includes all of the frequencies for the first SPW
        indices = [frequencies.index(h) for h in frequencies if (h < min_freq or h > max_freq)]
        number_of_freq += len(indices)

        if min(frequencies) < min_freq:
            min_freq = min(frequencies)
        if max(frequencies) > max_freq:
            max_freq = max(frequencies)

        sigma = 7 * RMS
        peaks, properties = find_peaks(kelvin_values, height=sigma, distance=9)

        #Excluding the peaks that have frequencies from overlapping SPWs
        kept_peaks = [frequencies[i] for i in peaks if i in indices]

        all_peaks.extend(kept_peaks)
        all_heights.extend([properties['peak_heights'][i] for i in range(len(properties['peak_heights'])) if i in indices])
        all_kelvin.extend([kelvin_values[i] for i in peaks if i in indices])
        all_rms.extend([RMS] * len(kept_peaks))

    snr = [(all_kelvin[h] / all_rms[h]) for h in range(len(all_rms))]
    GHz = number_of_freq / 1000
    peaks_per_GHz = round((len(all_peaks) / GHz), 5)
    avg = [peaks_per_GHz] * len(snr)
    final_peak_per_GHz.append(peaks_per_GHz)
    d = {'Frequencies (MHz)': all_peaks, 'Temperature (K)': all_kelvin, 'Signal to Noise Ratio': snr, 'Peak per GHz': avg}
    pd.DataFrame(data=d).to_csv(start_dir + 'SPARKS_spectra/Peak_Counts/' + str(round(sigma/RMS)) + 'Sigma/peaks_' + filename + '.dat', index=False)
    
data = {'Sources': all_file_names, 'Peaks per GHz': final_peak_per_GHz, 'Groups': group_numbers[::4]}
pd.DataFrame(data=data).to_csv(start_dir + 'SPARKS_spectra/' + str(round(sigma/RMS)) + '_sigma_peaks.dat', index=False)

print("--- Runtime: %s seconds ---" % round((time.time() - start_time), 2))