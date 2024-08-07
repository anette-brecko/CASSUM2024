import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import os
start_time = time.time()

source_info_dir = '/Users/elsabrecko/Desktop/CASSUM/SPARKS_spectra/source_info_NEW.txt'
start_dir = '/Users/elsabrecko/Desktop/CASSUM/SPARKS_spectra/GROUP'
methanol_df = pd.read_csv('/Users/elsabrecko/Desktop/CASSUM/SPARKS_spectra/Peak_Counts/methanol_lines_5Sigma.dat')
avg_offset = methanol_df.iloc[:, 9]
methanol_sources = methanol_df.iloc[:, 9]

#Define constants
df = pd.read_csv(source_info_dir)
k = 1.380649 * pow(10, -23)
c = 299792458
omega_conversion = math.pi / (4 * np.log(2))
arcsec_to_rad = 4.84814 * pow(10, -6)

#Pandas has seems to have difficulty reading in the source_info_NEW.txt file due to formatting, 
#so I'm using for loops and list comprehension to get the desired values

#I have 2 seperate lists for the protocluster/source names and omega_beam values

#corrected_frequency = frequency * (1 + vlsr / c)

split_list = [x for y in df.iloc[:, 0] for x in y.split()]
group_list = split_list[16::12]
name_list = []
omega_list = []
vlsr_list = []
cluster_count = 17
mm_count = 19
omega_count = 22
vlsr_count = 25
for i in range(len(group_list)):
    name = split_list[cluster_count] + split_list[cluster_count+1] + '-MM' + split_list[mm_count]
    name_list.append(name)
    theta_maj = float(split_list[omega_count]) * arcsec_to_rad
    theta_min = float(split_list[omega_count+1]) * arcsec_to_rad
    omega = theta_maj * theta_min * omega_conversion
    omega_list.append(omega)
    offset = 0
    if name in methanol_sources:
        offset = avg_offset[i]
    corrected_vlsr = float(split_list[vlsr_count]) + offset
    vlsr_list.append(corrected_vlsr * 1000)
    cluster_count += 12
    mm_count += 12 
    omega_count += 12
    vlsr_count += 12

#Finding central frequencies and computing Kelvin values from source files, which have the format XXXX_peak_continuum.line.dat
for i in range(len(group_list)):
    kelvin_dir = start_dir + str(group_list[i]) + '_K/'
    end_file = '_peak_continuum.line.dat'
    start_file = start_dir + str(group_list[i]) + '/' + str(name_list[i]) + '_'
    if os.path.exists(kelvin_dir):
        pass
    else:
        os.mkdir(kelvin_dir)
    for j in range(4):
        dir = start_file + str(j) + end_file
        dataframe = pd.read_csv(dir)
        split_list2 = [x for y in dataframe.iloc[:, 0] for x in y.split()]
        frequencies = list(map(float, split_list2[::2]))
        frequencies = [(i * 1000000) for i in frequencies]
        jy_beam = list(map(float, split_list2[1::2]))
        central_freq = (min(frequencies) + max(frequencies)) / 2
        conversion_factor = (c * c * pow(10, -26)) / (2 * omega_list[i] * k * central_freq * central_freq)
        
        #Multiply each Jy_beam value by conversion factor, change back to MHz, account for vLSR shift, write to new .dat file
        kelvin = [(i * conversion_factor) for i in jy_beam]
        frequencies = [(i / 1000000) for i in frequencies]
        frequencies = [(x * (1 + (vlsr_list[i] / c))) for x in frequencies]
        df2 = pd.DataFrame(data={'col0': frequencies, 'col1': kelvin})
        new_dir = kelvin_dir + str(name_list[i]) + '_' + str(j) + '_K.dat'
        df2.to_csv(new_dir, index=False)

print("--- Runtime: %s seconds ---" % round((time.time() - start_time), 5))