# theta_physiology
All code and statistics used in the paper titled "Theta mediated dynamics of human hippocampal-neocortical learning systems in memory formation and retrieval" published in Nature Communications in 2023. First author is Sandra Gattas.

Abstract:
Episodic memory is thought to arise as a function of dynamic interactions between the hippocampus and the neocortex, yet the mechanisms by which these interactions manifest have remained elusive. Here, using human intracranial recordings during a mnemonic discrimination task, we report that 4-5 Hz (overlapping with theta) power is differentially recruited during pattern separation compared to overgeneralization errors, and its phase supports interactions between the hippocampus and the neocortex when memories are being formed and correctly retrieved. Interactions were largely bidirectional, with small but significant net directional biases; a hippocampus-to-neocortex bias during the acquisition of new information that was subsequently correctly discriminated from similar stimuli, and a neocortex-to-hippocampus bias during the accurate discrimination of new stimuli from similar previously learned stimuli at retrieval. These results may be somewhat surprising given predictions from Complementary Learning Systems and similar models and provide new insights into the mechanisms by which the two learning systems (hippocampus and neocortex) may support episodic memory. The 4-5 Hz rhythm may facilitate the initial stages of information acquisition by neocortex during learning (i.e., hippocampus->neocortex information transfer) and the recall of stored information from cortex during retrieval (i.e., neocortex->hippocampus information transfer), in the context of a discrimination memory task. Future work should further probe these dynamics across different types of tasks and stimuli and computational models may need to be expanded accordingly to accommodate these findings.

In addition to the code base here, there are several other packages that are necessary to download:

MVGC_v1.0 (https://github.com/SacklerCentre/MVGC1) -  MVGC Multivariate Granger Causality toolbox
NLX2MAT (https://github.com/tjd2002/nlx2mat-win/) - Package to extract Neuralynx data to MATLAB
osort_v4.0 (https://github.com/shihchengyen/osort-v4-rel) - Spike detection and sorting MATLAB package
FieldTrip (https://github.com/fieldtrip/fieldtrip) - MATLAB Toolbox for MEG, EEG, and iEEG analysis (release 20181130)
Brainstorm 3 (https://github.com/brainstorm-tools/brainstorm3) - Brainstorm software for electrophysiology analysis.
