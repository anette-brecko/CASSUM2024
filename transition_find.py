from astroquery.linelists.cdms import CDMS
from astroquery.splatalogue import Splatalogue
import astropy.units as u
import pandas as pd
import pprint

start_dir = '/Users/elsabrecko/Desktop/CASSUM/SPARKS_spectra/GROUP'
df = pd.read_csv('/Users/elsabrecko/Desktop/CASSUM/source_info_converted.txt')
sources = df.iloc[:, 0]
groups = df.iloc[:, 5]
spw01_min = float('inf')
spw01_max = float('-inf')
spw23_min = float('inf')
spw23_max = float('-inf')

for i in range(len(sources)):
    file = start_dir + str(groups[i]) + '_K/' + sources[i] + '_K.dat'
    frequencies = pd.read_csv(file).iloc[:, 0]
    spw = int(file[len(file) - 7])
    current_min = min(frequencies)
    current_max = max(frequencies)
    
    if spw == 0 or spw == 1:
        if current_max > spw01_max:
            spw01_max = current_max
        elif current_min < spw01_min:
            spw01_min = current_min
            
    if spw == 2 or spw == 3:
        if current_max > spw23_max:
            spw23_max = current_max
        elif current_min < spw23_min:
            spw23_min = current_min

print('SPW01 Min: ' + str(spw01_min))
print('SPW01 Max: ' + str(spw01_max))
print('SPW23 Min: ' + str(spw23_min))
print('SPW23 Max: ' + str(spw23_max))

columns = ('chemical_name',
               'orderedfreq',
               'aij',
               'upper_state_energy_K')
spw01 = Splatalogue.query_lines(min_frequency=(spw01_min * u.MHz), max_frequency=(spw01_max * u.MHz), chemical_name = ' SO ',
                                energy_max=1000, energy_type='eu_k', line_lists=['CDMS'])
mask01 = spw01['aij'] > -8.0
spw01 = spw01[columns]
methanol_freq_01 = spw01['orderedfreq']
spw01.pprint(max_width = 300)
#pd.DataFrame(data=spw01.to_pandas()).to_csv('/Users/elsabrecko/Desktop/CASSUM/spw01_methyl_formate.dat', index=False)

spw23 = Splatalogue.query_lines(min_frequency=(spw23_min * u.MHz), max_frequency=(spw23_max * u.MHz), chemical_name = ' SO ',
                                energy_max=1000, energy_type='eu_k', line_lists=['CDMS'])
mask = spw23['aij'] > -8.0
spw23 = spw23[columns]
methanol_freq_23 = spw23['orderedfreq']
spw23.pprint(max_width = 300)
#pd.DataFrame(data=spw23.to_pandas()).to_csv('/Users/elsabrecko/Desktop/CASSUM/spw23_methyl_formate.dat', index=False)