clear;
clc;
clf;

AllChinIDs= 369; %[366 369]; %[369 369 370 358 360 366 367];

% tStart = .033; tEnd = .167;
% tStart = .235; tEnd = .325;
% tStart = .235; tEnd = .58;
% tStart = .180; tEnd = .720;
% tStart = .38; tEnd = .48;
% tStart = .48; tEnd = .58;
% tStart = .578; tEnd = .668;
% tStart = .74; tEnd = .85;
% tStart = .89; tEnd = 1.19;
% tStart = .92; tEnd = 1.02;
% tStart = 1.02; tEnd = 1.15;
tStart = 0; tEnd = 1.3;


saveFigs= 1;
rmsENV= zeros(length(AllChinIDs), 5); % where columns are [S | SN | N]
rmsTFS= zeros(length(AllChinIDs), 5); % where columns are [S | SN | N]


for chinVar= 1:length(AllChinIDs)
    chinID= AllChinIDs(chinVar);
    RootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
    
    allFiles= dir([RootDataDir '*' num2str(chinID) '*SFR*']);
    
    if isempty(allFiles)
        error('No dir. what to do?');
    elseif length(allFiles)>1
        fprintf('there are multiple dirs. \n');
        
        for dirVar= 1:length(allFiles)
            fprintf('(%d)-%s\n', dirVar, allFiles(dirVar).name);
        end
        chosen_dir_num= input('Which one? \n');
        
        data_dir= [RootDataDir allFiles(chosen_dir_num).name filesep];
    else
        data_dir= [RootDataDir allFiles.name filesep];
    end
    
    fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/Segment_%.0fto%.0fms%s', tStart*1e3, tEnd*1e3, filesep);
    if ~isdir(fig_save_dir)
        mkdir(fig_save_dir);
    end
    
    s_files= dir([data_dir 'a*_S_*']);
    sn_files= dir([data_dir 'a*_SN_*']);
    n_files= dir([data_dir 'a*_N_*']);
    
    %% clean speech
    s_data_cell= cell(length(s_files), 2);
    nPairs_actual= nan(length(s_files), 1);
    for sfile_var=1:length(s_files)
        temp_data= load([data_dir s_files(sfile_var).name]);
        temp_data = temp_data.data;
        s_data_cell{sfile_var, 1}= temp_data.AD_Data.AD_Avg_PO_V{1};
        s_data_cell{sfile_var, 2}= temp_data.AD_Data.AD_Avg_NP_V{1};
        
        nPairs_actual(sfile_var)= temp_data.Stimuli.RunLevels_params.nPairs_actual;
    end
    
    s_atten=temp_data.Stimuli.atten_dB;
    
    [s_sig, fs_sig]= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_09_11-Q362_AN_PTS/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');
    t_sig= (1:length(s_sig))/fs_sig;
    
    s_data_pos= zeros(1, length(s_data_cell{sfile_var,1}));
    s_data_neg= zeros(1, length(s_data_cell{sfile_var,2}));
    fs_data= temp_data.Stimuli.RPsamprate_Hz;
    
    
    for i=1:length(s_files)
        s_data_pos= s_data_pos + s_data_cell{i, 1}*nPairs_actual(i)/sum(nPairs_actual);
        s_data_neg= s_data_neg + s_data_cell{i, 2}*nPairs_actual(i)/sum(nPairs_actual);
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
    
    
    %% plot Clean speech
    figHan= 1;
    fName_S= sprintf('Q%d_SFR_speech_alone_atten%.0fdB', chinID, s_atten);
    ttlStr= sprintf('clean speech, %d reps, atten%.0fdB', sum(nPairs_actual), s_atten);
    create_panel_plot_include_hilbert(figHan, fs_sig, s_sig, fs_data, s_data_pos_filt, s_data_neg_filt, s_data_env, s_data_tfs, tStart, tEnd, fig_save_dir, fName_S, ttlStr, saveFigs);
    
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
    
    %% plot noisy speech
    figHan= 2;
    fName_SN= sprintf('Q%d_SFR_noisy_speech_0dB_SNR_atten%.0fdB', chinID, sn_atten);
    ttlStr= sprintf('noisy speech, %d reps, atten%.0fdB', sum(sn_nreps), sn_atten);
    create_panel_plot_include_hilbert(figHan, fs_sig, sn_sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, sn_data_env, sn_data_tfs, tStart, tEnd, fig_save_dir, fName_SN, ttlStr, saveFigs);
    
    
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
    
    %% plot
    
    %% plot Clean speech
    figHan= 3;
    fName_N= sprintf('Q%d_SFR_noise_alone_atten%.0fdB', chinID, n_atten);
    ttlStr= sprintf('noise-alone, %d reps, atten%.0fdB', sum(n_nreps), n_atten);
    create_panel_plot_include_hilbert(figHan, fs_sig, n_sig, fs_data, n_data_pos_filt, n_data_neg_filt, n_data_env, n_data_tfs, tStart, tEnd, fig_save_dir, fName_N, ttlStr, saveFigs);
    
    %% other comps
    rmsENV(chinVar, :)= [rms(s_data_env), rms(sn_data_env) rms(n_data_env) s_atten contains(data_dir, 'PTS')];
    rmsTFS(chinVar, :)= [rms(s_data_tfs) rms(sn_data_tfs) rms(n_data_tfs) s_atten contains(data_dir, 'PTS')];
end

rmsRatio= rmsENV;
rmsRatio(:, 1:3)= rmsENV(:, 1:3) ./ rmsTFS(:, 1:3);

%%
function curFilt= get_filter(fs_data)
N_bp_half= 4;
HalfPowerFrequency1=0.5;
HalfPowerFrequency2=4e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end