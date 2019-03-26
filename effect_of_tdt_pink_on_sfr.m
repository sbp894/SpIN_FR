% tests effect of tdt-generated pink noise on SFR

clear;
clc;

chinID= 370;
RootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
fig_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/tdt_pink_SFR/';
stimDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/stim/';
saveFigs= 1;

if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end

data_dir= dir([RootDataDir '*Q' num2str(chinID) '*SFR_pink*']);
data_dir= [RootDataDir data_dir.name filesep];

allfiles= dir([data_dir 'a*SFR*.mat']);
allfiles= allfiles(~(contains({allfiles.name}, 'latency') | contains({allfiles.name}, 'artifact')));

all_snrs= cell2mat(cellfun(@(x) str2double(strrep(x(regexp(x, 'snr_')+4 : regexp(x, '_atn')-1), 'm', '-')), {allfiles.name}, 'uniformoutput', false));

tStart= 0; tEnd= 1.3;
tStart_whole= 0; tEnd_whole= 1.3;
fig_save_dir_subdir= [fig_save_dir sprintf('t%.0fto%.0f_ms/', tStart*1e3, tEnd*1e3)];


raw_power_env= nan(length(all_snrs), 1);
frac_power_env= nan(length(all_snrs), 1);
raw_power_tfs= nan(length(all_snrs), 1);
frac_power_tfs= nan(length(all_snrs), 1);

stim_fName= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2017_11_02-Q325_AN_NH/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav';
[sig, fs_sig]= audioread(stim_fName);

snrVar= 1;
[s_data_pos_filt, s_data_neg_filt]= get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs);



for snrVar= 2:length(all_snrs)
    curSNR= all_snrs(snrVar);
    
    [sn_data_pos_filt, sn_data_neg_filt, fs_data]= get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs);
    
    figHan= 1;
    fName= strrep(sprintf('Q%d_nh_SNR%d_sn', chinID, curSNR), '-', 'm');
    ttlStr= sprintf('Q%d,NH,SNR %d', chinID, curSNR);
    
    PSD_struct= create_panel_plot_s_vs_sn(figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, s_data_pos_filt, s_data_neg_filt, tStart, tEnd, ...
        fig_save_dir_subdir, fName, ttlStr, saveFigs);
    
    
    raw_power_env(snrVar)= PSD_struct.raw.ENV.SN;
    raw_power_tfs(snrVar)= PSD_struct.raw.TFS.SN;
    frac_power_env(snrVar)= PSD_struct.frac.ENV.SN;
    frac_power_tfs(snrVar)= PSD_struct.frac.TFS.SN;
    
    if snrVar==length(all_snrs)
        raw_power_env(1)= PSD_struct.raw.ENV.S;
        raw_power_tfs(1)= PSD_struct.raw.TFS.S;
        frac_power_env(1)= PSD_struct.frac.ENV.S;
        frac_power_tfs(1)= PSD_struct.frac.TFS.S;
    end
    
end

%%
figure(11);
clf;
mrkSize= 20;
lw= 4;
fSize= 20;
[~, sort_inds]= sort(all_snrs);

xlabel_str= cellfun(@(x) num2str(x), num2cell([all_snrs(sort_inds)]), 'uniformoutput', false);


subplot(211);
hold on;
plot(1:size(raw_power_env,1), raw_power_env(sort_inds), 'v', 'markersize', mrkSize, 'linew', lw);
plot(1:size(raw_power_tfs,1), raw_power_tfs(sort_inds), '^', 'markersize', mrkSize, 'linew', lw);
% legend('ENV', 'TFS');
grid on;
ylabel('$RAW_{power} (dB)$', 'interpreter', 'latex');
set(gca,'fontsize', fSize, 'xtick', 1:size(raw_power_env,1), 'xticklabel', xlabel_str);
title(['NH- Q' num2str(chinID)]);

subplot(212);
hold on;
plot(1:size(frac_power_env,1), frac_power_env(sort_inds), 'v', 'markersize', mrkSize, 'linew', lw);
plot(1:size(frac_power_tfs,1), frac_power_tfs(sort_inds), '^', 'markersize', mrkSize, 'linew', lw);
legend('ENV', 'TFS', 'location', 'southeast');
grid on;
ylabel('$FRAC_{power} (dB)$', 'interpreter', 'latex');
set(gca,'fontsize', fSize, 'xtick', 1:size(raw_power_env,1), 'xticklabel', xlabel_str);
xlabel('HP Masking SNR (dB)');
set(gcf, 'units', 'inches', 'position', [1 1 10 6]);

fName_summary= sprintf('Q%d_nh_pink_summary', chinID);
saveas(gcf, [fig_save_dir_subdir fName_summary], 'png');


function curFilt= get_filter(fs_data)
N_bp_half= 10;
HalfPowerFrequency1=50;
HalfPowerFrequency2=0.5e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end

function [sn_data_pos_filt, sn_data_neg_filt, fs_data]= get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs)
curSNR= all_snrs(snrVar);
curFile= allfiles(snrVar).name;

fprintf('Using %s for %d dB SNR\n', curFile, curSNR);


temp= load([data_dir curFile]);

sn_data_pos= temp.data.AD_Data.AD_Avg_PO_V{1};
sn_data_neg= temp.data.AD_Data.AD_Avg_NP_V{1};
fs_data= temp.data.Stimuli.RPsamprate_Hz;

initialRampDur= 50e-3;
ramp_nSamples= round(initialRampDur*fs_data);
bpFilt= get_filter(fs_data);

sn_rampHamming= hamming(2*ramp_nSamples)';
sn_rampVector= [sn_rampHamming(1:ramp_nSamples), ones(1, length(sn_data_pos)-length(sn_rampHamming)) sn_rampHamming(ramp_nSamples+1:end)];

sn_data_pos_filt= filtfilt(bpFilt, sn_data_pos.*sn_rampVector);
sn_data_neg_filt= filtfilt(bpFilt, sn_data_neg.*sn_rampVector);
end