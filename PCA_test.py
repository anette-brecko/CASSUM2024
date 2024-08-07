import numpy as np
import pandas as pd
import glob
from sklearn.decomposition import PCA
from sklearn.metrics import root_mean_squared_error
from sklearn.utils.multiclass import type_of_target
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astroquery.linelists.cdms import CDMS
from astroquery.splatalogue import Splatalogue
import astropy.units as u
from scipy.signal import find_peaks


##############
### STEP 1 ###
##############
# Load all synthetic spectra (= training set) and pre-process them
# All spectra need to have same length (i.e. same frequency grid)
# All spectra need to be re-aligned to account for frequency/velocity shift

spectra = []
SPW = '1'
mypath = "/Users/elsabrecko/Desktop/CASSUM/PCA/"   
source_info = pd.read_csv('/Users/elsabrecko/Desktop/CASSUM/source_info_converted.txt')
file_list = glob.glob(mypath + 'SPW' + SPW + '/*.dat')
methanol_freq = [335133.570, 335582.017, 336438.224, 336865.149, 337135.853, 345903.916, 345919.260, 346202.719]

# Determine a common frequency grid
common_frequency = None

# To plot all synthetc spectra
plt.figure(figsize=(10, 6)) 

for file_name in file_list:
    # Read the file line by line, skipping the first line if it starts with "!"
    with open(file_name, 'r') as file:
        lines = file.readlines()
        if lines[0].startswith('!'):
            lines = lines[1:]
        data = np.loadtxt(lines)
        
    frequency = data[:, 0]
    intensity = data[:, 1]
    
    #Finding peaks in synthetic spectrum, using arbitrarily low height value
    peaks, properties = find_peaks(intensity, height=0.5, distance=9)
    peak_freq = [frequency[i] for i in peaks]
    
    #Getting methanol transitions for the given frequency range of the SPW
    columns = ('chemical_name',
               'orderedfreq',
               'aij',
               'upper_state_energy_K')
    query = Splatalogue.query_lines(min_frequency=(min(frequency) * u.MHz), max_frequency=(max(frequency) * u.MHz), chemical_name = ' CH3OH ',
                                energy_max=1000, energy_type='eu_k', line_lists=['CDMS'])
    mask = query['aij'] > -8.0
    methanol_freq = query['orderedfreq']
    for i in methanol_freq:
        if any((i - 50) <= x < (i + 50) for x in peak_freq):
            val = list(j for j in peak_freq if (i - 50) <= j <= (i + 50))
            offset = i - val[0]
            frequency = [(f + offset) for f in frequency]
            break
    
    if common_frequency is None:
        common_frequency = frequency
    else:
        min_freq = min(min(common_frequency), min(frequency))
        max_freq = max(max(common_frequency), max(frequency))
        common_frequency = np.linspace(min_freq, max_freq, num=len(common_frequency))
  
    # Interpolate to the common frequency grid
    interp_func = interp1d(frequency, intensity, kind='linear', fill_value="extrapolate")
    intensity_interp = interp_func(common_frequency)


    # Plot each spectrum
    plt.plot(common_frequency, intensity_interp, label=file_name)

    spectra.append(intensity_interp)
    

# Convert to numpy array
spectra = np.array(spectra)

# Verify the shape of the spectra array
print(f'Spectra shape: {spectra.shape}')

# Plot
plt.xlabel('Frequency (MHz)')
plt.ylabel('Intensity (K)')
plt.title('Synthetic Spectra')
plt.grid(True)
plt.tight_layout()
plt.show()

# Here it shows that (as you know already) there can be velocity/frequency offsets between the synthetic spectra
# Therefore, you will need to re-align the synthetic spectra within this step one




##############
### STEP 2 ###
##############
# Apply PCA to the training set

n_components = 5  # Number of principal components to keep, this can be edited, try different values
pca = PCA(n_components=n_components)
pca.fit(spectra)

eigenvectors = pca.components_ 
eigenvalues = pca.explained_variance_

print(f'Eigenvectors:\n{eigenvectors}') 
print(f'Eigenvalues:\n{eigenvalues}')

# Transform the spectra using PCA
spectra_pca = pca.transform(spectra)

# Plot the explained variance
plt.figure(figsize=(8, 6))
plt.plot(np.cumsum(pca.explained_variance_ratio_), marker='o')
plt.xlabel('Number of Components')
plt.ylabel('Cumulative Explained Variance') #cumulative proportion of total variance explained by the components up to that point
plt.title('Explained Variance by PCA Components') #also called "Scree Plot"
plt.grid()
plt.show()

# The scree plot helps decide how many principal components to retain for dimensionality reduction while preserving most of the variance.
# e.g. if 95% of the variance is explained by the first 3 components, you might choose to use just those 3 components instead of the entire original set of variables, i.e. reducing the dimensionality of the data from the original number of features to just three components.
# Retaining too many components can lead to overfitting



##############
### STEP 3 ###
##############
# Test on observed spectrum 

# Load observed spectrum
#Find group number of the source, access the Kelvin file for this source and SPW
source = '351.1542+0.7073-MM23_' + SPW
group = str(source_info.loc[source_info.iloc[:, 0] == source, 'group number'].values[0])

#Getting the 
path_obs = "/Users/elsabrecko/Desktop/CASSUM/SPARKS_Spectra/GROUP" + group + "_K/" + source + "_K.dat"
new_spectrum_data = pd.read_csv(path_obs)
frequency_obs = list(new_spectrum_data.iloc[:, 0])
intensity_obs = list(new_spectrum_data.iloc[:, 1])

# Interpolate observed spectrum to the common frequency grid
interp_func = interp1d(frequency_obs, intensity_obs, kind='linear', fill_value="extrapolate")
intensity_obs = interp_func(common_frequency)

# Convert to numpy array
spectra_obs = np.array([intensity_obs])

# Verify the shape of the spectra array
print(f'Spectra_obs shape: {spectra_obs.shape}')

# Project the observed spectrum onto the PCA components
new_spectrum_pca = pca.transform(intensity_obs.reshape(1, -1))

# Print projected spectrum (PCA scores) 
print(f'Projected observed spectrum (PCA scores):\n{new_spectrum_pca}')

# Reconstruct the spectrum from the principal components
reconstructed_spectrum = pca.inverse_transform(new_spectrum_pca)  # = approximation of the original spectrum created using a reduced set of principal components
for i in range(len(reconstructed_spectrum[0])):
    curr = reconstructed_spectrum[0][i]

# Plot the original and reconstructed spectrum for comparison
plt.figure(figsize=(10, 6))
plt.plot(common_frequency, intensity_obs, label='Original Observed Spectrum')
plt.plot(common_frequency, reconstructed_spectrum.flatten(), label='Reconstructed Spectrum')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Intensity (K)')
plt.legend()
plt.title('Original vs Reconstructed Spectrum')
plt.show()



##############
### STEP 4 ###
##############
# Determine whether methanol lines are detected
# To detect specific lines, compare the reconstructed spectrum to the original one, only where the reconstructed spectrum intensity is not zero
non_zero_indices = reconstructed_spectrum.flatten() != 0

# Define a threshold for line detection
threshold = 0.9 # Edit threshold to decide which percentage of the line intensity you want 

# Compute difference between observed spectrum and reconstructed spectrum, only for channels where reconstructed spectrum is != 0
detected_lines = np.abs(intensity_obs[non_zero_indices] / reconstructed_spectrum.flatten()[non_zero_indices]) > threshold
line_frequencies = common_frequency[non_zero_indices][detected_lines]

#Print all channels for which > XX% of the expected methanol lines are found in the observed spectrum
# where XX is determined by the selected threshold
#print("Detected line frequencies (MHz):", line_frequencies)
count = 0
average_ratio = 0
RMS = source_info.loc[source_info.iloc[:, 0] == source, 'rms (Kelvin)'].values[0]

for i in line_frequencies:
    if min(frequency_obs) <= i <= max(frequency_obs):
        curr = round(i, 2)
        rounded = [round(j, 2) for j in frequency_obs]
        val = list(h for h in rounded if (curr - 0.5) <= h <= (curr + 0.5))
        index = rounded.index(val[0])
        observed_intensity = intensity_obs[index]
        reconstructed_intensity = reconstructed_spectrum[0][index]
        if reconstructed_intensity >= (5*RMS) and observed_intensity >= (5*RMS):
            minimum = min(reconstructed_intensity, observed_intensity)
            maximum = max(reconstructed_intensity, observed_intensity)
            count += 1
            ratio = minimum / maximum
            average_ratio += ratio
            if ratio >= 0.9:
                print("Ratio: " + str(ratio))

print("Final number of frequencies over 5 sigma (for both observed + reconstructed): " + str(count))
if count == 0:
    print("No methanol peaks detected")
else:
    rmse = root_mean_squared_error(y_true=intensity_obs[non_zero_indices], y_pred=reconstructed_spectrum.flatten()[non_zero_indices])
    print("Root Mean Squared Error: " + str(rmse))
    print("Percent Error: " + str(rmse / len(non_zero_indices)))
    average_ratio = average_ratio / count
    print("Average accuracy: " + str(average_ratio))