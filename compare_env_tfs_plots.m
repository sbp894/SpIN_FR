clear;

s369 = load('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Artifact_Removed_FFR/SP-2018_11_05-Q369_SFR_NH_pilot2/a0002_FFR_SNRenvSSN_Stim_S_P_nType0_atten10.mat');
s370= load('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Artifact_Removed_FFR/SP-2018_11_05-Q370_SFR_NH/a0002_FFR_SNRenvSSN_Stim_S_P_nType0_atten10.mat');
s373= load('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Artifact_Removed_FFR/SP-2019_07_17-Q373_SFR_NH/a0008_FFR_SNRenvSSN_Stim_S_P_atn10.mat');

s369_pos= s369.data.AD_Data.AD_Avg_PO_V{1};
s370_pos= s370.data.AD_Data.AD_Avg_PO_V{1};
s373_pos= s373.data.AD_Data.AD_Avg_PO_V{1};

s369_neg= s369.data.AD_Data.AD_Avg_NP_V{1};
s370_neg= s370.data.AD_Data.AD_Avg_NP_V{1};
s373_neg= s373.data.AD_Data.AD_Avg_NP_V{1};

fs = s369.data.Stimuli.RPsamprate_Hz;
t= (1:length(s369_neg))/fs;

bfFilt= designfilt('bandpassiir','FilterOrder',4, ...
'HalfPowerFrequency1',100,'HalfPowerFrequency2',1500, ...
'SampleRate', fs);


lw= 2;
figure(1); clf;
ax(1)= gca;
hold on
plot(t, filtfilt(bfFilt, s369_pos+s369_neg), 'LineWidth', lw)
plot(t, filtfilt(bfFilt, s370_pos+s370_neg), 'LineWidth', lw)
plot(t, filtfilt(bfFilt, s373_pos+s373_neg), 'LineWidth', lw)
title('ENV')
grid on

figure(2); clf;
hold on
ax(2)= gca;
plot(t, filtfilt(bfFilt, s369_pos-s369_neg), 'LineWidth', lw)
plot(t, filtfilt(bfFilt, s370_pos-s370_neg), 'LineWidth', lw)
plot(t, filtfilt(bfFilt, s373_pos-s373_neg), 'LineWidth', lw)
title('TFS')
grid on

linkaxes(ax)
lg= legend('369 ', '373', '370');
lg.FontSize= 20;

