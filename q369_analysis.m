clear;
clc;
clf;

data_dir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_10_30-Q369_SFR_pilot1/';
fig_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Pilot Figures/';
lw= 2;

s_files= dir([data_dir 'a*_S_*']);
sn_files= dir([data_dir 'a*_SN_*']);
n_files= dir([data_dir 'a*_N_*']);

%% clean speech
sn_data_cell= cell(length(s_files), 2);
nPairs_actual= nan(length(s_files), 1);
for sfile_var=1:length(s_files)
    temp_data= load([data_dir s_files(sfile_var).name]);
    temp_data = temp_data.data;
    sn_data_cell{sfile_var, 1}= temp_data.AD_Data.AD_Avg_PO_V{1};
    sn_data_cell{sfile_var, 2}= temp_data.AD_Data.AD_Avg_NP_V{1};
    
    nPairs_actual(sfile_var)= temp_data.Stimuli.RunLevels_params.nPairs_actual;
end

[s_sig, fs_sig]= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_09_11-Q362_AN_PTS/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');
t_sig= (1:length(s_sig))/fs_sig;

s_data_pos= zeros(1, length(sn_data_cell{sfile_var,1}));
s_data_neg= zeros(1, length(sn_data_cell{sfile_var,2}));
fs_data= temp_data.Stimuli.RPsamprate_Hz;
t_data= (1:length(s_data_pos))/fs_data;


for i=1:length(s_files)
    s_data_pos= s_data_pos + sn_data_cell{sfile_var, 1}*nPairs_actual(sfile_var)/sum(nPairs_actual);
    s_data_neg= s_data_neg + sn_data_cell{sfile_var, 2}*nPairs_actual(sfile_var)/sum(nPairs_actual);
end


N_bp_half= 4;
HalfPowerFrequency1=0.5;
HalfPowerFrequency2=1e3;
initialRampDur= 50e-3;
ramp_nSamples= round(initialRampDur*fs_data);
rampHamming= hamming(2*ramp_nSamples)';
rampVector= [rampHamming(1:ramp_nSamples), ones(1, length(s_data_pos)-length(rampHamming)) rampHamming(ramp_nSamples+1:end)];
curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
s_data_pos_filt= filtfilt(curFilt, s_data_pos.*rampVector).*rampVector;
s_data_neg_filt= filtfilt(curFilt, s_data_neg.*rampVector).*rampVector;
s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;

%% plot
figure(1);
clf;
Ashift= 0.125;
fSize= 16;
ax(1)= plot(t_data, Ashift + s_data_pos_filt, 'linew', lw);
hold on;
ax(2)= plot(t_data, Ashift + s_data_neg_filt, 'linew', lw);
plot(nan, nan);
ax(3)= plot(t_data, s_data_env, 'linew', lw);
ax(3)= plot(t_data,- Ashift + s_data_tfs, 'linew', lw);
ax(4)= plot(t_sig, - 2*Ashift + s_sig*0.2, 'k', 'linew', lw);
grid on;
legend(ax, '+ve pol', '-ve pol', 'sum of pols', 'stim');
set(gca, 'fontsize', fSize);
xlabel('time (sec)');
ylabel('Amp');
title(sprintf('clean speech, 65 dB SPL, %d reps', sum(nPairs_actual)));


set(gcf, 'units', 'normalized', 'position', [0 0 1 1]);
fName_S= 'Q369_SFR_speech_alone';
saveas(gcf, [fig_save_dir fName_S]);
saveas(gcf, [fig_save_dir fName_S], 'tiff');

%% noisy speech
temp_data= load([data_dir sn_files.name]);
temp_data = temp_data.data;
sn_data_pos= temp_data.AD_Data.AD_Avg_PO_V{1};
sn_data_neg= temp_data.AD_Data.AD_Avg_NP_V{1};
sn_nreps= temp_data.Stimuli.RunLevels_params.nPairs_actual;

sn_sig= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_09_11-Q362_AN_PTS/Signals/MH/SNRenv/SNR_0/SSN_Stim0dB_SN_P.wav');

fs_data= temp_data.Stimuli.RPsamprate_Hz;
t_data= (1:length(s_data_pos))/fs_data;


sn_data_pos_filt= filtfilt(curFilt, sn_data_pos.*rampVector).*rampVector;
sn_data_neg_filt= filtfilt(curFilt, sn_data_neg.*rampVector).*rampVector;
sn_data_env= (sn_data_pos_filt+sn_data_neg_filt)/2;
sn_data_tfs= (sn_data_pos_filt-sn_data_neg_filt)/2;

%% plot
figure(2);
clf;
Ashift= 0.05;
fSize= 16;
ax(1)= plot(t_data, Ashift + sn_data_pos_filt, 'linew', lw);
hold on;
ax(2)= plot(t_data, Ashift + sn_data_neg_filt, 'linew', lw);
plot(nan, nan);
ax(3)= plot(t_data, sn_data_env, 'linew', lw);
ax(3)= plot(t_data, - Ashift  + sn_data_tfs, 'linew', lw);
ax(4)= plot(t_sig, - 2*Ashift + sn_sig*Ashift, 'k', 'linew', lw);
grid on;
legend(ax, '+ve pol', '-ve pol', 'sum of pols', 'stim');
set(gca, 'fontsize', fSize);
xlabel('time (sec)');
ylabel('Amp');
title(sprintf('noisy-speech, 0 dB SNR, 65 dB SPL, %d reps', sn_nreps));

set(gcf, 'units', 'normalized', 'position', [0 0 1 1]);
fName_SN= 'Q369_SFR_noisy_speech';
saveas(gcf, [fig_save_dir fName_SN]);
saveas(gcf, [fig_save_dir fName_SN], 'tiff');

%% noisy speech
temp_data= load([data_dir n_files.name]);
temp_data = temp_data.data;
n_data_pos= temp_data.AD_Data.AD_Avg_PO_V{1};
n_data_neg= temp_data.AD_Data.AD_Avg_NP_V{1};
n_nreps= temp_data.Stimuli.RunLevels_params.nPairs_actual;

n_sig= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_09_11-Q362_AN_PTS/Signals/MH/SNRenv/SNR_0/SSN_Stim0dB_N_P.wav');

fs_data= temp_data.Stimuli.RPsamprate_Hz;
t_data= (1:length(s_data_pos))/fs_data;


n_data_pos_filt= filtfilt(curFilt, n_data_pos.*rampVector).*rampVector;
n_data_neg_filt= filtfilt(curFilt, n_data_neg.*rampVector).*rampVector;
n_data_env= (n_data_pos_filt+n_data_neg_filt)/2;
n_data_tfs= (n_data_pos_filt-n_data_neg_filt)/2;

%% plot
figure(3);
clf;
Ashift= 0.1;
fSize= 16;
ax(1)= plot(t_data, Ashift + n_data_pos_filt, 'linew', lw);
hold on;
ax(2)= plot(t_data, Ashift + n_data_neg_filt, 'linew', lw);
plot(nan, nan);
ax(3)= plot(t_data, n_data_env, 'linew', lw);
ax(3)= plot(t_data,- Ashift + n_data_tfs, 'linew', lw);
ax(4)= plot(t_sig, - 2*Ashift + n_sig*Ashift, 'k', 'linew', lw);
grid on;
legend(ax, '+ve pol', '-ve pol', 'sum of pols', 'stim');
set(gca, 'fontsize', fSize);
xlabel('time (sec)');
ylabel('Amp');
title(sprintf('noise alone, 65 dB SPL, %d reps', sn_nreps));

set(gcf, 'units', 'normalized', 'position', [0 0 1 1]);
fName_N= 'Q369_SFR_noise_alone';
saveas(gcf, [fig_save_dir fName_N]);
saveas(gcf, [fig_save_dir fName_N], 'tiff');