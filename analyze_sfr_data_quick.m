% run('/media/parida/DATAPART1/Matlab/ExpData/NelData/SP-2017_05_15-Q313_SFR_Pilot/a0007_FFR_SNRenvStim0dB_SN_P_nType0_atten25.m');
% run('/media/parida/DATAPART1/Matlab/ExpData/NelData/SP-2017_05_15-Q314_SFR_Pilot/a0003_FFR_SNRenvStim0dB_SN_P_nType1_atten25.m');
run('/media/parida/DATAPART1/Matlab/ExpData/NelData/SP-2017_05_15-Q313_SFR_Pilot/a0007_FFR_SNRenvStim6dB_SN_P_nType0_atten25.m')
data= ans;
[stim, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/stimSetStationary/Stim6dB_SN_P.wav');
t_stim= (1:length(stim))/fs;

np= data.AD_Data.AD_Avg_NP_V{1};
po= data.AD_Data.AD_Avg_PO_V{1};
t_efr= (1:length(po))/data.Stimuli.RPsamprate_Hz;

lw=2;
clf
plot(t_efr, po, 'linew', lw)
hold on
plot(t_efr, np, 'linew', lw)
plot(t_stim, -.02+ .01*stim, 'k')
grid on
ylim([-.03 .02]);
xlim([0 2])