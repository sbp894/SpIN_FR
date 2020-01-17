clear;

load('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2019_05_08-Q371_SFRpink500Hz_NH/p0002_SFR_pink_S_snr_120_atn10.mat')
all_data=data.AD_Data.AD_All_V';


pos_data = cell2mat(all_data(1:2:end));
t_data= (1:size(pos_data, 2)) / data.Stimuli.RPsamprate_Hz;

N_bp_half= 10;
HalfPowerFrequency1=50;
HalfPowerFrequency2=5e3;
curFilt= designfilt('bandpassiir','FilterOrder',4, ...
'HalfPowerFrequency1',50,'HalfPowerFrequency2',3e3, ...
'SampleRate',data.Stimuli.RPsamprate_Hz);


nf_pos= reshape(pos_data(randsample(numel(pos_data), numel(pos_data))), size(pos_data,1), size(pos_data,2));

pos_data_filt= filtfilt(curFilt, mean(pos_data, 1));
pos_nf_filt= filtfilt(curFilt, mean(nf_pos, 1));
pos_nf_filt_sub= filtfilt(curFilt, mean(pos_data(1:2:end, :)-pos_data(2:2:end, :), 1));

figure(1);
clf;
hold on;
plot(t_data, pos_data_filt, 'linew', 2)
plot(t_data, pos_nf_filt, 'linew', 1.5);
plot(t_data, pos_nf_filt_sub, 'linew', 1.5);
ylim([-.1 .1])


figure(2);
clf;
fs_data= data.Stimuli.RPsamprate_Hz;
plot_dpss_psd(pos_data_filt, fs_data, 'nw', 5)
hold on
plot_dpss_psd(pos_nf_filt, fs_data, 'nw', 5)
plot_dpss_psd(pos_nf_filt_sub, fs_data, 'nw', 5)
xlim([50 3e3]);


durStim= 1.3; 
indData= t_data<=durStim;
indNF= t_data>durStim;

% figure(3);
% clf;
% plot_dpss_psd(pos_data_filt(indData), fs_data, 'nw', 5)
% hold on
% plot_dpss_psd(pos_nf_filt(indNF), fs_data, 'nw', 5)
% xlim([50 3e3]);
% 
% figure(1)
% plot(t_data(indData), pos_data_filt(indData), '--')
% plot(t_data(indNF), pos_data_filt(indNF), '--')
% 
