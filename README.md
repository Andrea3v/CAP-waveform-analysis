# CAP-waveform-analysis
analysis of partial CAP recorded waveforms, returning area, peak amplitude, and latency
The Script works on .mat exported data from Patchmaster v>=2.9x
data can be exported from both channels (V-mon on channel 2). 
if only one channel, *A/D is found, the value of pre-amplification used is required. 

the script uses .1 - .9 ms as a range for waveform offset & automatically detects the range for integration of the CAP area by: 
1) calculating the mean of the baseline waveforms
2a) using the prompt 'end of stim' [ms] to detect the beginning of the evoked CAP on the average waveforms
2b) If the rise of the evoked CAP response cannot be found, the user is required to add it manually by inserting a data tip to determine the left range for integration
3) using the peak latency of the baseline multiplied by a factor, e.g. 125%, to set the right end of the integration 

the range for calc of peak amplitude and latency has to be inserted manually by adding two data tips that determine the range where to find the waveform MAX

the script will output the analyzed data as .xlsx in the subfolder 'Analysis' along with the final figure and the variables used for analysis
The script won't run on exported data which might contain traces with both, one AND two channels, exporting only one type of data at the time
