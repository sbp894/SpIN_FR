% For pilot, gonna record 250 reps. Approximately good PSD estimate. (SP)
% 11/7/18

% clear;
clf;
clc;

% s_data= load('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_10_31-Q358_PTS_SFR/p0002_FFR_SNRenvSSN_Stim_S_P_nType0_atten10.mat');
s_data= load('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_10_31-Q360_PTS_SFR/p0002_FFR_SNRenvSSN_Stim_S_P_nType0_atten10.mat');

s_data= s_data.data;
pos_s_data= (s_data.AD_Data.AD_All_V(1:2:end))';
neg_s_data= (s_data.AD_Data.AD_All_V(2:2:end))';
fs_data= s_data.Stimuli.RPsamprate_Hz;

onsetZero_time= 2e-3;
onsetZero_inds= 1:round(onsetZero_time*fs_data);

all_nReps= [500 350 250 150];

N_bp_half= 4;
HalfPowerFrequency1=0.5;
HalfPowerFrequency2=1e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);

all_temp_data_pos= cell(length(all_nReps), 1);
all_temp_data_neg= cell(length(all_nReps), 1);

ax= nan(length(all_nReps), 1);
bx= nan(length(all_nReps), 1);
for nRepVar= 1:length(all_nReps)
    nRep= all_nReps(nRepVar);
    
    temp_data_pos= mean(cell2mat(pos_s_data(1:nRep)), 1);
    temp_data_pos= filtfilt(curFilt, temp_data_pos);
    temp_data_pos(onsetZero_inds)= 0;
    all_temp_data_pos{nRepVar}= temp_data_pos;
    
    temp_data_neg= mean(cell2mat(neg_s_data(1:nRep)), 1);
    temp_data_neg= filtfilt(curFilt, temp_data_neg);
    temp_data_neg(onsetZero_inds)= 0;
    all_temp_data_neg{nRepVar}= temp_data_neg;
    
    tData= (1:length(temp_data_pos))/fs_data;
    
    figure(1);
    ax(nRepVar)= subplot(length(all_nReps), 2, 2*nRepVar-1);
    hold on;
    plot(tData, temp_data_pos, '-', tData, temp_data_neg, '-');
    pos_residue= all_temp_data_pos{1} - temp_data_pos;
    neg_residue= all_temp_data_neg{1} - temp_data_neg;
    plot(tData, pos_residue, '-.k', tData, neg_residue, '--g');
    title(sprintf('Error = %.1f', 100*( rms(pos_residue)/rms(all_temp_data_pos{1}) + rms(neg_residue)/rms(all_temp_data_neg{1}))/2));
    ylabel(sprintf('nReps=%d', nRep));
    
    bx(nRepVar)= subplot(length(all_nReps), 2, 2*nRepVar);
    hold on;
    plot_dpss_psd(temp_data_pos, fs_data, 'NW', 4);
    plot_dpss_psd(pos_residue, fs_data, 'NW', 4);
end
xlabel('time (sec)');
linkaxes(ax);
linkaxes(bx);