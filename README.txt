README
Estimating stress drop using the Generalized Inversion Technique (GIT)

Step 1
INPUT: all files in RC_beta directory
Script: stations_step1.py
Purpose: create list in event folder of all stations which we have records from for that event
OUTPUT: list of stations (stations.csv) in event directory

Step 2
INPUT: all files in RC_beta directory
Script: stations_step2.py
Purpose: determine which of the stations exists in *.txt files containing site response information
OUTPUT: list of stations and whether or not (boolean) we have site response info from them (stations_inv.csv) in event directory

Step 3
INPUT: all files in RC_beta directory
Script: instrument_corrections.py
Purpose: take data and detrend and apply instrument correction to it.  The *.xml files in the data set do not contain response information for all stations, and so a list of uncorrected files is created, and these are not saved in resulting data.  
OUTPUT: all files in corrected directory

Step 4
INPUT: all files in corrected directory
Script: compute_spectra_meters.py
Purpose: take data and compute power spectra binned by frequency
OUTPUT: all files in record_spectra directory

Step 5
INPUT: all files in record_spectra directory
Script: secondo_meters.py
Purpose: take record spectra and discard all files from stations with less than 3(?) records, and run Andrews inversion to obtain event and station spectra.  
OUTPUT: all files in Andrews_inversion directory

Step 6
INPUT: all files in Andrews_inversion directory
Script: findBrune_trapezoids.py
Purpose: find amplitude constraint in form of event spectra that looks most Brune like
OUTPUT: constraint function file in constraint directory


Step 7
INPUT: all files in Andrews_inversion directory
Script: secondo_constrain.py
Purpose: apply constraint function to all event and station spectra
OUTPUT: constrained spectra in Andrews_inversion_constrained directory

Step 8
INPUT: all files in Andrews_inversion_constrained directory
Script: fitBrune.py
Purpose: take constrained event spectra, preform non-linear least squares inversion and find a best fit moment and corner frequency for the brune model for each event, then use those parameters to calculate stress drop for each event
OUTPUT: dataframe containing for each event: event id, magnitude, catalog moment, guessed fc, various goodness-of-fit parameters, best fit fc, best fit moment, and stress drop

