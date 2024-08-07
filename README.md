# CASSUM 2024
Code used for my project in the CASSUM summer program, focusing on automating the analysis 
of radio astrononmical spectra for molecular indentification and characterization.

[Program information](https://cosmicorigins.space/cassum-vico24)<br />
[Final presentation](https://drive.google.com/file/d/1K9DS--Lo1gwMuQSdogHVzZdRzedzqb2x/view)

This project used data from the SPARKS observing program (Search for high-mass protostars with ALMA revealed up 
to kiloparsec scales), from which the individual data files are not yet publicly available.

### File Structure:
1. `Jybeam_to_Kelvin.py`<br />
Intensities of the source data files were in units of Jansky per beam, which is dependent on
the radio telescope from which the data is gathered. This script read in all of the source
files and converted the intensities to units of Kelvin.
2. `find_peaks.ipynb`<br />
Reads in provided information about the sources through source_info_converted.txt, which includes
RMS, rest velocity of the sources, size of the radio telescope beam, etc. and finds spectral peaks
in the data above a particular threshold (i.e. 3σ, 5σ). Plots peaks in all four spectral
windows and displays these results using matplotlib.
3. `find_all_peaks.py`<br />
Finds molecular transition peaks in all of the ~200 source files, and saves the frequencies and
intensities of the peaks in each source as a .csv file.
4. `transition_find.py`<br />
Finds the overall minimum and maximum frequencies across all sources, and queries
the Splatalogue database for astronomical spectroscopy for particular molecules
to see which transitions of a given molecule fall within this frequency range. Saves
the resulting frequencies, Einstein coefficients, and upper state energy to a .dat file.
5. `PCA_test.py`<br />
Uses a small number of sythentic spectra generated via LTE fitting particular molecule
and the PCA package in scikit-learn in order to make a reconstructed approximation of
the observed data. The accuracy of this reconstruction allows us to determine whether or
not a transition for that particular molecule is detected within the observed spectra.
