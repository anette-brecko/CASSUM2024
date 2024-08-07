import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import pandas as pd

#Want to keep track of the frequencies that have been recorded in a previous iteration of the loop

start_dir = '/Users/elsabrecko/Desktop/CASSUM/'
meth_peaks = pd.read_csv(start_dir + 'SPARKS_Spectra/Peak_Counts/methanol_5Sigma_5MHz.dat')
end_file = '_K.dat'

def plot_spectra(filename):
    source_info = pd.read_csv(start_dir + 'source_info_converted.txt')
    col = source_info.iloc[:, 0]
    #all_kelvin has all temperature values for the peaks
    all_kelvin = []
    #all_heights has the heights of the peaks in the order of which they're found
    all_heights = []
    #all_peaks has the frequencies of all the peaks
    all_peaks = []
    #There's one RMS value per SPW
    all_rms = []
    min_freq = float('inf')
    max_freq = float('-inf')
    avg_peaks = 0
    #number_of_freq keeps track of how many frequencies are covered (minus overlap) in order to calcualte average peaks/GHz
    number_of_freq = 0
    spw0 = [347478.251, 347493.965, 347628.340, 347617.010, 347604.659, 348065.967, 348049.886]
    spw1 = [345974.664, 345985.381, 346001.616, 346675.644, 346687.469]
    spw2 = [336028.165, 336111.324, 336351.390, 336918.184, 336889.198]
    spw3 = [334031.781]
    methyl_formate = spw0 + spw1 + spw2 + spw3
    methanol = [335133.57, 335582.017, 336438.224, 336865.149, 337135.853, 345903.916, 345919.26, 346202.719]
    
    
    for i in range(4):
        name = filename + '_' + str(i)
        RMS = source_info.loc[col == name, 'rms (Kelvin)'].values[0]
        
        #Will eventually need to modify this so the group number is read in from source_info_converted.txt
        group = source_info.loc[source_info.iloc[:, 0] == name, 'group number'].values[0]
        dir = start_dir + 'SPARKS_spectra/GROUP' + str(group) + '_K/' + name + end_file
        df = pd.read_csv(dir)
        frequencies = list(df.iloc[:, 0])
        kelvin_values = list(df.iloc[:, 1])
        
        #Get index of frequencies that have not yet been covered by a different SPW
        #Should include all of the frequencies for the first SPW
        indices = [frequencies.index(h) for h in frequencies if (h < min_freq or h > max_freq)]
        number_of_freq += len(indices)
        
        if min(frequencies) < min_freq:
            min_freq = min(frequencies)
        if max(frequencies) > max_freq:
            max_freq = max(frequencies)

        sigma = 5 * RMS
        peaks, properties = find_peaks(kelvin_values, height=sigma, distance=9)
        fig, ax = plt.subplots(figsize=(20, 5))
        
        #Excluding the peaks that have frequencies from overlapping SPWs
        kept_peaks = [frequencies[i] for i in peaks if i in indices]
        
        all_peaks.extend(kept_peaks)
        all_heights.extend([properties['peak_heights'][i] for i in range(len(properties['peak_heights'])) if i in indices])
        all_kelvin.extend([kelvin_values[i] for i in peaks if i in indices])
        all_rms.extend([RMS] * len(kept_peaks))

        # Change plot color to black
        ax.step(frequencies, kelvin_values, color='black')
        ax.scatter(np.array(frequencies)[peaks], np.array(kelvin_values)[peaks], color='red')
        ax.axhline(y = sigma, color= 'b', linestyle='--')
        for j in methyl_formate:
            for q in kept_peaks:
                if min(frequencies) <= j <= max(frequencies):
                    if (j - 1.5) <= q <= (j + 1.5):
                        ax.axvline(x = j, color= 'green', linestyle='--')
        for j in range(len(methanol)):
            curr = methanol[j]
            for q in kept_peaks:
                if min(frequencies) <= curr <= max(frequencies):
                    if (curr - 5) <= q <= (curr + 5):
                        ax.axvline(x = curr, color= 'blue', linestyle='--')
        
        ax.set_xlabel('Frequency (MHz)')
        ax.set_ylabel('Temperature (Kelvin)')

        # Calculate peaks per GHz (convert MHz to GHz first)
        freq_range = (max(frequencies) / 1000) - (min(frequencies) / 1000)
        peaks_per_GHz = round((len(kept_peaks) / freq_range), 5)
        avg_peaks += peaks_per_GHz

        # Display number of peaks and peaks per GHz on the plot
        ax.text(0.02, 0.95, f'Number of Peaks: {len(peaks)}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.4, 0.95, f'{name}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.9, 0.95, f'Threshold: {round(sigma/RMS)}σ', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.02, 0.90, f'Peaks per GHz: {peaks_per_GHz}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        plt.tight_layout()
        plt.show()
    
    snr = [(all_kelvin[i] / all_rms[i]) for i in range(len(all_rms))]
    GHz = number_of_freq / 1000
    avg = [round((len(all_peaks) / GHz), 5)] * len(snr)
    d = {'Frequencies (MHz)': all_peaks, 'Temperature (K)': all_kelvin, 'Signal to Noise Ratio': snr, 'Peak per GHz': avg}
    #pd.DataFrame(data=d).to_csv(start_dir + 'peaks_' + filename + '.dat', index=False)

#input_filemame should just be the group, source, and MM number
#The code loops through the individual spectral window files for each source
input_filename = '329.1835-0.3147-MM3'
plot_spectra(input_filename)