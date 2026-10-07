READ ME:

This data set is associated with the following manuscript published in PNAS (2016):

Title:
A face-selective ventral occipito-temporal map of the human brain with intracerebral potentials.

Authors: 
Jacques Jonas, Corentin Jacques, Joan Liu-Shuang, Hélène Brissart, Sophie Colnat-Coulbois, Louis Maillard and Bruno Rossion.



--- 1. Electrophysiological data ---

Each folder contains the data from one participant. The data consist in recordings of intracerebral EEG in epileptic patient viewing visual stimuli. Each channel corresponds to the recording (i.e. voltage over time) of one intracerebral contact and additional physiological markers (e.g. ECG). 

In each folder there is data from two conditions (see file name): 
1. face_periodic: the main condition of the paper: periodic presentation of faces
2. face_nonperiodic: control condition: non periodic presentation of faces

Each condition consist of two files (.mat and .lw5) that can be open in Letswave 5, a free software for EEG data analysis (http://www.nocions.org/letswave/) running on Matlab. 

Prior to this step, sEEG data (TRC files, Micromed) have been pre-processed in the following way:

1. Imported in Letswave 5, which has been used for the subsequent steps
2. Low-pass filtered with a cut-off at 30Hz using a butterworth low-pass filter (Order 4)
3. Each recording sequence was segmented from 2 seconds after the onset of the sequence until ~65 seconds, so that the number of bins in the sequence allows an exact integer number of cycles of the 1.2 Hz face frequency. 
4. Sequences for each condition were merged in the same file: 2 or 4 sequences for the 'face_periodic' and 1 or 2 for the 'face_non_periodic'. 


Further data transform and analyses (FFT, etc.) have not been performed.


Letswave files can be open in Matlab without 'Letswave'. They consist of one header file (the .lw5 file) and one data file (the .mat file). 
The dimensions of the 'data' variable are organized as follow:

-dim 1: sequences/epochs
-dim 2: channels
-dim 3: unused
-dim 4: unused
-dim 5: unused
-dim 6: time


--- 2. Recording contacts coordinates ---

In addition, the locations of each intracerebral channel both in the Talairach and the MNI space are reported in a txt file in each participant's folder. Note that these coordinates were used only for display purposes in the current paper. Anatomical definitions were performed on the individual participant's anatomy. 



--- 3. Set of stimuli used ---


The full set of face and objects used in the study are included in a separate folder (/Stimuli/).











