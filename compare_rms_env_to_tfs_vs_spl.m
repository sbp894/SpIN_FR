clear;
clc;
clf;

AllChinIDs= 368; %[369 369 370 358 360 366 367];

tStart = 0; tEnd = 1.3;

fig_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/RMS_ratio_env_to_tfs/';
if ~isdir(fig_save_dir)
    mkdir(fig_save_dir);
end


NumOfWindows= 100; windDur= 50e-3;
% NumOfWindows= 75; windDur= 100e-3;
windTimes_data= nan(NumOfWindows, 2);
windTimes_data(:,1)= (tEnd-windDur)*rand(NumOfWindows, 1);
windTimes_data(:,2)= windTimes_data(:,1) + windDur;

tFinal= 1.5;
windTimes_nf= nan(NumOfWindows, 2);
windTimes_nf(:,1)= tEnd + 5e-5 + (tFinal-tEnd)*rand(NumOfWindows, 1);
windTimes_nf(:,2)= windTimes_nf(:,1) + windDur;


saveFigs= 1;
rmsENV= zeros(length(AllChinIDs), 6); % where columns are [chinID | S | SN | N | atten | hearing]
rmsTFS= zeros(length(AllChinIDs), 6); % where columns are [chinID | S | SN | N | atten | hearing]
rmsRatio_speech= nan(length(AllChinIDs), NumOfWindows);
rmsRatio_noisefloor= nan(length(AllChinIDs), NumOfWindows);
segSPLs_speech= nan(length(AllChinIDs), NumOfWindows);

for chinVar= 1:length(AllChinIDs)
    chinID= AllChinIDs(chinVar);
    RootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
    
    allFiles= dir([RootDataDir '*' num2str(chinID) '*SFR*']);
    allFiles= allFiles(~contains({allFiles.name}, 'pink'));
    
    if isempty(allFiles)
        error('No dir. what to do?');
    elseif length(allFiles)>1
        warning('there are multiple dirs. choosing last one');
        warning('bad!!! ');
        if chinVar==1
            data_dir= [RootDataDir allFiles(1).name filesep];
        else
            data_dir= [RootDataDir allFiles(2).name filesep];
        end
    else
        data_dir= [RootDataDir allFiles.name filesep];
    end
    
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
    
    s_atten=temp_data.Stimuli.atten_dB;
    
    [s_sig, fs_sig]= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_09_11-Q362_AN_PTS/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');
    t_sig= (1:length(s_sig))/fs_sig;
    
    s_data_pos= zeros(1, length(sn_data_cell{sfile_var,1}));
    s_data_neg= zeros(1, length(sn_data_cell{sfile_var,2}));
    fs_data= temp_data.Stimuli.RPsamprate_Hz;
    
    
    for i=1:length(s_files)
        s_data_pos= s_data_pos + sn_data_cell{sfile_var, 1}*nPairs_actual(sfile_var)/sum(nPairs_actual);
        s_data_neg= s_data_neg + sn_data_cell{sfile_var, 2}*nPairs_actual(sfile_var)/sum(nPairs_actual);
    end
    
    initialRampDur= 50e-3;
    ramp_nSamples= round(initialRampDur*fs_data);
    rampHamming= hamming(2*ramp_nSamples)';
    rampVector= [rampHamming(1:ramp_nSamples), ones(1, length(s_data_pos)-length(rampHamming)) rampHamming(ramp_nSamples+1:end)];
    curFilt= get_filter(fs_data);
    
    s_data_pos_filt= filtfilt(curFilt, s_data_pos.*rampVector).*rampVector;
    s_data_neg_filt= filtfilt(curFilt, s_data_neg.*rampVector).*rampVector;
    s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
    s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;
    
    
    %% noisy speech
    temp_data= load([data_dir sn_files.name]);
    temp_data = temp_data.data;
    sn_data_pos= temp_data.AD_Data.AD_Avg_PO_V{1};
    sn_data_neg= temp_data.AD_Data.AD_Avg_NP_V{1};
    sn_nreps= temp_data.Stimuli.RunLevels_params.nPairs_actual;
    sn_sig= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_09_11-Q362_AN_PTS/Signals/MH/SNRenv/SNR_0/SSN_Stim0dB_SN_P.wav');
    
    fs_data= temp_data.Stimuli.RPsamprate_Hz;
    sn_atten= temp_data.Stimuli.atten_dB;
    
    sn_data_pos_filt= filtfilt(curFilt, sn_data_pos.*rampVector).*rampVector;
    sn_data_neg_filt= filtfilt(curFilt, sn_data_neg.*rampVector).*rampVector;
    sn_data_env= (sn_data_pos_filt+sn_data_neg_filt)/2;
    sn_data_tfs= (sn_data_pos_filt-sn_data_neg_filt)/2;
    
    
    %% noisy speech
    temp_data= load([data_dir n_files.name]);
    temp_data = temp_data.data;
    n_data_pos= temp_data.AD_Data.AD_Avg_PO_V{1};
    n_data_neg= temp_data.AD_Data.AD_Avg_NP_V{1};
    n_nreps= temp_data.Stimuli.RunLevels_params.nPairs_actual;
    
    n_sig= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_09_11-Q362_AN_PTS/Signals/MH/SNRenv/SNR_0/SSN_Stim0dB_N_P.wav');
    
    n_atten= temp_data.Stimuli.atten_dB;
    fs_data= temp_data.Stimuli.RPsamprate_Hz;
    t_data= (1:length(s_data_pos))/fs_data;
    
    
    n_data_pos_filt= filtfilt(curFilt, n_data_pos.*rampVector).*rampVector;
    n_data_neg_filt= filtfilt(curFilt, n_data_neg.*rampVector).*rampVector;
    n_data_env= (n_data_pos_filt+n_data_neg_filt)/2;
    n_data_tfs= (n_data_pos_filt-n_data_neg_filt)/2;
    
    
    for windVar=1:NumOfWindows
        curInds2use_data= t_data>windTimes_data(windVar, 1) & t_data<windTimes_data(windVar, 2);
        rmsRatio_speech(chinVar, windVar)= rms(s_data_env(curInds2use_data))/rms(s_data_tfs(curInds2use_data));
        
        curInds2use_nf= t_data>windTimes_nf(windVar, 1) & t_data<windTimes_nf(windVar, 2);
        rmsRatio_noisefloor(chinVar, windVar)= rms(s_data_env(curInds2use_nf))/rms(s_data_tfs(curInds2use_nf));
        
        inds2use_stim= t_sig>windTimes_data(windVar, 1) & t_sig<windTimes_data(windVar, 2);
        segSPLs_speech(chinVar, windVar)= 20*log10(rms(s_sig(inds2use_stim))/(20e-6));
    end
    
    %% other comps
    rmsENV(chinVar, :)= [chinID rms(s_data_env), rms(sn_data_env) rms(n_data_env) s_atten contains(data_dir, 'PTS')];
    rmsTFS(chinVar, :)= [chinID rms(s_data_tfs) rms(sn_data_tfs) rms(n_data_tfs) s_atten contains(data_dir, 'PTS')];
end

rmsRatio= rmsENV;
rmsRatio(:, 1:3)= rmsENV(:, 1:3) ./ rmsTFS(:, 1:3);
segSPLs_speech= unique(segSPLs_speech, 'rows');

rmsRatio_noisefloor= rmsRatio_noisefloor(:);
spl_noisefloor_rand= min(segSPLs_speech) + range(segSPLs_speech)*rand(size(rmsRatio_noisefloor));

%%
figure(1);
clf;
hold on;
mrkSize= 14;
fSize= 16;
lw=2;
ax(1)= plot(nan, nan, 'bd', 'markersize', mrkSize, 'linew', lw);
ax(2)= plot(nan, nan, 'kd', 'markersize', mrkSize, 'linew', lw);
ax(3)= plot(nan, nan, 'ro', 'markersize', mrkSize, 'linew', lw);
% ax(4)= plot(nan, nan, 'g.', 'markersize', mrkSize, 'linew', lw);
% ax(2)= plot(nan, nan, 'mo');

% plot(spl_noisefloor_rand, rmsRatio_noisefloor, '.g', 'markersize', mrkSize, 'linew', lw);
for chinVar= 1:length(AllChinIDs)
    if rmsENV(chinVar, end)==0 % NH
        if rmsENV(chinVar, end-1)==10 % 10 dB Atten
            plot(segSPLs_speech, rmsRatio_speech(chinVar, :), 'bd', 'markersize', mrkSize, 'linew', lw);
        elseif rmsENV(chinVar, end-1)==25 % 25 dB Atten
            plot(segSPLs_speech, rmsRatio_speech(chinVar, :), 'kd', 'markersize', mrkSize, 'linew', lw);
        end
        
    elseif rmsENV(chinVar, end)==1 % PTS
        if rmsENV(chinVar, end-1)==10 % 10 dB Atten
            plot(segSPLs_speech, rmsRatio_speech(chinVar, :), 'ro', 'markersize', mrkSize, 'linew', lw);
        elseif rmsENV(chinVar, end-1)==25 % 25 dB Atten
            plot(segSPLs_speech, rmsRatio_speech(chinVar, :), 'mo', 'markersize', mrkSize, 'linew', lw);
        end
    end
end
grid on;
xlabel('Segment dB SPL');
ylabel('rms(ENV/TFS)');
title(sprintf('Window = %.0f ms, #Windows= %.0f', windDur*1e3, NumOfWindows));
set(gca, 'fontsize', fSize);
legend(ax, 'NH (10 attn)', 'NH (25 attn)', 'PTS (10 attn)');
set(gcf, 'units', 'normalized', 'position', [0 0 1 1]);
axis tight;

fName= sprintf('NH_vs_PTS_rms_ratio_of_env_tfs_%.0fms_win', windDur*1e3);
saveas(gcf, [fig_save_dir fName], 'tiff');

function curFilt= get_filter(fs_data)
N_bp_half= 4;
HalfPowerFrequency1=0.5;
HalfPowerFrequency2=1e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end