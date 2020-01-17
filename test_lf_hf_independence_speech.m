% Testing if LF and HF in voiced segments are independent or not
clear;
clf;
clc;



%% define directory with wav-files 


% First start with one wav-file 
% read speech
[sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');

% test for voicing: extract voiced segments 


% estimate HF and LF power in some time bin 


% plot correlations for voiced segments: see if correlations are different
% for unvoiced segments 