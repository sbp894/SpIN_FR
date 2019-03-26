% Script to create FFR panel plots per chin
clear;
clc;

AllChinIDs= [358 360 366 367 369 370]; % [366 369];

[sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');
t_sig= (1:length(sig))/fs;
stim_dur= length(sig)/fs;

restricted_time= [0 .72; 0.88 1.2];
% restricted_time= [0 1.3];
saveAgain= 1;

windowLength= 64e-3;
fracOverLap= .75;
NW= 1.5; % NW= dur * f_res => f_res= NW/dur. If NW=1.5, dur=50ms, f_res= 30 Hz

fracSlide= 1-fracOverLap;
tSlide= fracSlide*windowLength;
nSegs= 1 + floor((stim_dur-windowLength)/(tSlide));

dpss_yRange= 100;
plot_dpss_each_iter= false;
lw= 2;
lw2= 3;
nSProws=2;
fSize= 16;
nSPcols=1;

t_latency= 5e-3;

ratio_hf_to_lf_audio= nan(nSegs, 1);
audio_freq_band_low= [50 500];
audio_freq_band_high= [500 3000];

data_f0_related_band= [50 500];
ratio_f0_env_to_tfs= nan(nSegs, 1);

for chinVar= 1:length(AllChinIDs)
    restart_fig= 1;
    figure(2);
    clf;
    co= get(gca, 'colororder');
    
    if plot_dpss_each_iter
        figure(1);
        clf;
        hold on;
    end
    
    chinID= AllChinIDs(chinVar);
    
    if chinID==191
        warning('Don''t have data in response to speech, Collected data to N twice. Recollect data');
        pause
    end
    
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
        if contains(allFiles(chosen_dir_num).name, 'NH')
%             chosen_dir_num= chosen_dir_num;
            postFix= 'NH';
        elseif contains(allFiles(chosen_dir_num).name, 'PTS')
%             chosen_dir_num= 1;
            postFix= 'PTS';
%         elseif contains(allFiles(chosen_dir_num).name, 'NH')
% %             chosen_dir_num= 2;
%             postFix= 'NH';
        end
        fprintf('choosing %d\n', chosen_dir_num);
        
        
        data_dir= [RootDataDir allFiles(chosen_dir_num).name filesep];
    else
        data_dir= [RootDataDir allFiles.name filesep];
        if contains(data_dir, 'NH')
            postFix= 'NH';
        elseif contains(data_dir, 'PTS')
            postFix= 'PTS';
        end
    end
    
    fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/moving_segment_analysis_vowel/');
    if ~isdir(fig_save_dir)
        mkdir(fig_save_dir);
    end
    
    s_files= dir([data_dir 'a*_S_*']);
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
    
    initialRampDur= 20e-3;
    ramp_nSamples= round(initialRampDur*fs_data);
    rampHamming= hamming(2*ramp_nSamples)';
    rampVector= [rampHamming(1:ramp_nSamples), ones(1, length(s_data_pos)-length(rampHamming)) rampHamming(ramp_nSamples+1:end)];
    curFilt= get_filter(fs_data);
    
    s_data_pos_filt= filtfilt(curFilt, s_data_pos.*rampVector).*rampVector;
    s_data_neg_filt= filtfilt(curFilt, s_data_neg.*rampVector).*rampVector;
    s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
    s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;
    t_data= (1:length(s_data_env))/fs_data;
    
    A_stim= max([s_data_env s_data_tfs]);
    A_shift= 1.2*A_stim;
    
    for segVar= 1:nSegs
        seg_t_start=  (segVar-1)*tSlide;
        seg_ind_start= max(1, round(seg_t_start*fs));
        seg_t_end= seg_t_start + windowLength;
        seg_ind_end= min(length(sig), round(seg_t_end*fs));
        
        seg_t_mid= (seg_t_start+seg_t_end)/2;
        seg_in_rest_time= (seg_t_mid>restricted_time(:,1)) & (seg_t_mid<restricted_time(:,2));
        if any(seg_in_rest_time)
            
            
            seg_stim= sig(seg_ind_start:seg_ind_end);
            seg_t= t_sig(seg_ind_start:seg_ind_end);
            
            cur_data_inds= t_data>(seg_t_start+t_latency) & t_data<(seg_t_end+t_latency);
            cur_data_env= s_data_env(cur_data_inds);
            cur_data_tfs= s_data_tfs(cur_data_inds);
            cur_t_data= t_data(cur_data_inds);
            
            if plot_dpss_each_iter
                subplot(nSProws, nSPcols, 1);
                if restart_fig
                    restart_fig= 0;
                    clf;
                    subplot(nSProws, nSPcols, 1);
                    hold on;
                    plot(t_sig, -A_shift + A_stim*sig, 'color', co(1,:)); % very first time
                    plot(t_data, A_shift + s_data_env, 'color', co(4,:)); % very first time
                    plot(t_data, s_data_tfs, 'color', co(5,:)); % very first time
                    
                elseif exist('bx', 'var')
                    for bxVar=1:length(bx)
                        delete(bx(bxVar));
                    end
                else
                    hold on;
                    plot(t_sig, -A_shift + A_stim*sig, 'color', co(1,:)); % very first time
                    plot(t_data, A_shift + s_data_env, 'color', co(4,:)); % very first time
                    plot(t_data, s_data_tfs, 'color', co(5,:)); % very first time
                end
                
                bx(1)= plot(seg_t, -A_shift + A_stim*seg_stim, 'r', 'linew', lw);
                bx(2)= plot(cur_t_data, A_shift + cur_data_env, 'r', 'linew', lw);
                bx(3)= plot(cur_t_data, cur_data_tfs, 'r', 'linew', lw);
                
                if exist('h2', 'var')
                    delete(h2);
                end
                h2= subplot(nSProws, nSPcols, 2);
                yyaxis left;
            end
            
            [Pxx_sig_dB, freq_stim, ~]= plot_dpss_psd(seg_stim, fs, 'NW', NW, 'plot', plot_dpss_each_iter);
            Pxx_sig= 10.^(Pxx_sig_dB/10);
            if plot_dpss_each_iter
                ylim([max(Pxx_sig_dB)+10-dpss_yRange max(Pxx_sig)+10])
            end
            
            sig_freq_inds_low= freq_stim>audio_freq_band_low(1) & freq_stim<audio_freq_band_low(2);
            sig_freq_inds_high= freq_stim>audio_freq_band_high(1) & freq_stim<audio_freq_band_high(2);
            Pxx_sig_low= sum(Pxx_sig(sig_freq_inds_low));
            Pxx_sig_high= sum(Pxx_sig(sig_freq_inds_high));
            ratio_hf_to_lf_audio(segVar)= (Pxx_sig_high/range(audio_freq_band_high))/(Pxx_sig_low/range(audio_freq_band_low));
            
            if plot_dpss_each_iter
                yyaxis right;
                hold on;
            end
            [Pxx_env_dB, ~, px_env]= plot_dpss_psd(cur_data_env, fs_data, 'NW', NW, 'plot', plot_dpss_each_iter); %#ok<*ASGLU>
            Pxx_env= 10.^(Pxx_env_dB/10);
            
            [Pxx_tfs_dB, freq_data, px_tfs]= plot_dpss_psd(cur_data_tfs, fs_data, 'NW', NW, 'plot', plot_dpss_each_iter);
            Pxx_tfs= 10.^(Pxx_tfs_dB/10);
            if plot_dpss_each_iter
                set(px_env, 'color', co(4,:), 'linestyle', '-');
                set(px_tfs, 'color', co(5,:), 'linestyle', '-');
            end
            
            data_freq_inds= freq_data>data_f0_related_band(1) & freq_data<data_f0_related_band(2);
            Pxx_env_seg= sum(Pxx_env(data_freq_inds));
            Pxx_tfs_seg= sum(Pxx_tfs(data_freq_inds));
            
            ratio_f0_env_to_tfs(segVar)= Pxx_env_seg/Pxx_tfs_seg;
            
            if plot_dpss_each_iter
                ylim1= max([Pxx_env_dB;Pxx_tfs_dB])+10-dpss_yRange;
                ylim2= ylim1 + dpss_yRange;
                ylim([ylim1 ylim2])
            end
            %         pause(.1);
        else
            %             ratio_f0_env_to_tfs(segVar)= nan;
        end
    end
    %%
    ratio_hf_to_lf_audio_est= sort(ratio_hf_to_lf_audio);
    mdl= fitlm(ratio_hf_to_lf_audio, ratio_f0_env_to_tfs);
    c_m = mdl.Coefficients.Estimate;
    ratio_f0_env_to_tfs_est= c_m(1)+ c_m(2)*ratio_hf_to_lf_audio_est;
    
    [ratio_f0_env_to_tfs_max_vals, ratio_f0_env_to_tfs_max_inds]= sort(ratio_f0_env_to_tfs, 'descend','MissingPlacement','last');
    nMax2Plot= 1;
    ratio_f0_env_to_tfs_max_inds= ratio_f0_env_to_tfs_max_inds(1:nMax2Plot); % Plot only nMax2Plot segs
    seg_ind_start= nan(nMax2Plot, 1);
    seg_ind_end= nan(nMax2Plot, 1);
    for maxVar=1:length(ratio_f0_env_to_tfs_max_inds)
        seg_t_start=  (ratio_f0_env_to_tfs_max_inds(maxVar)-1)*tSlide;
        seg_ind_start(maxVar)= max(1, round(seg_t_start*fs));
        seg_t_end= seg_t_start + windowLength;
        seg_ind_end(maxVar)= min(length(sig), round(seg_t_end*fs));
    end
    
    figure(2);
    subplot(211);
    hold on
    plot(t_sig, -A_shift + A_stim*sig, 'color', co(1,:)); % very first time
    for maxVar=1:nMax2Plot
        plot(t_sig(seg_ind_start(maxVar):seg_ind_end(maxVar)), -A_shift + A_stim*sig(seg_ind_start(maxVar):seg_ind_end(maxVar)), 'r'); % very first time
        plot( [t_sig(seg_ind_start(maxVar)) t_sig(seg_ind_start(maxVar))], [-2*A_shift  2*A_shift], 'r', 'linew', lw);
        plot( [t_sig(seg_ind_end(maxVar)) t_sig(seg_ind_end(maxVar))], [-2*A_shift  2*A_shift], 'r', 'linew', lw);
        text(t_sig(round((seg_ind_start(maxVar)+seg_ind_end(maxVar))/2)), -1.1*A_shift + A_stim*sig(round((seg_ind_start(maxVar)+seg_ind_end(maxVar))/2)), sprintf('%d', maxVar));
    end
    plot(t_data, A_shift + s_data_env, 'color', co(4,:)); % very first time
    plot(t_data, s_data_tfs, 'color', co(5,:)); % very first time
    title(sprintf('Q%d, atten=%.0f dB, %s', chinID, s_atten, postFix));
    
    lg_hans(1)= plot(nan, nan, 'color', co(1,:));
    lg_hans(2)= plot(nan, nan, 'color', co(4,:));
    lg_hans(3)= plot(nan, nan, 'color', co(5,:));
    lg= legend(lg_hans, 'sig', 'env', 'tfs');
    grid on
    
    figName= sprintf('Q%d_atten%.0fdB_%s', chinID, s_atten, postFix);
    set(gca, 'fontsize', fSize);
    xlabel('time(sec)');
    set(lg, 'fontsize', 12, 'location', 'southeast');
    axis tight;
    
    subplot(223);
    yyaxis left;
    [Pxx_sig, freq_stim, px_stim]= plot_dpss_psd(sig, fs, 'NW', NW);
    ylim([max(Pxx_sig)+10-dpss_yRange max(Pxx_sig)+10])
    
    yyaxis right;
    hold on;
    [Pxx_env, ~, px_env]= plot_dpss_psd(s_data_env, fs_data, 'NW', NW);
    set(px_env, 'color', co(4,:), 'linestyle', '-');
    [Pxx_tfs, freq_data, px_tfs]= plot_dpss_psd(s_data_tfs, fs_data, 'NW', NW);
    ylim([max([Pxx_env;Pxx_tfs])+10-dpss_yRange max([Pxx_env;Pxx_tfs])+10])
    set(px_tfs, 'color', co(5,:), 'linestyle', '-');
    lg= legend('sig', 'env', 'tfs');
    set(lg, 'fontsize', 12, 'location', 'southwest');
    set(gca, 'fontsize', fSize, 'ytick', []);
    ylabel('');
    title('');
    xlabel('Freq (Hz)');
    
    subplot(224);
    mrkSize= 16;
    hold on;
    plot(ratio_hf_to_lf_audio, ratio_f0_env_to_tfs, 'o', 'markersize', mrkSize, 'linew', lw);
    for maxVar=1:nMax2Plot
        plot(ratio_hf_to_lf_audio(ratio_f0_env_to_tfs_max_inds(maxVar)), ratio_f0_env_to_tfs(ratio_f0_env_to_tfs_max_inds(maxVar)), 'ro', 'markersize', mrkSize, 'linew', lw2);
        text(ratio_hf_to_lf_audio(ratio_f0_env_to_tfs_max_inds(maxVar)), ratio_f0_env_to_tfs(ratio_f0_env_to_tfs_max_inds(maxVar)), sprintf('%d', maxVar));
    end
    plot(ratio_hf_to_lf_audio_est, ratio_f0_env_to_tfs_est, '-', 'linew', lw);
    ytick_vals= [.1 1 10 100];
    set(gca, 'xscale', 'log', 'yscale', 'log', 'ytick', ytick_vals);
    ylim([min(ytick_vals) max(ytick_vals)]);
    
%     xlabel(sprintf('Audio ratio (norm) HF (%.1f-%.1f kHz) to LF (%.2f-%.2f kHz)', audio_freq_band_high(1)/1e3, audio_freq_band_high(2)/1e3, audio_freq_band_low(1)/1e3, audio_freq_band_low(2)/1e3));
%     ylabel(sprintf('ENV/TFS in FFR in %.2f-%.2f kHz', data_f0_related_band(1)/1e3, data_f0_related_band(2)/1e3));
    
    xlabel(sprintf('$Carrier( ^{MF=[.5-3kHz]}/_{LF[50-500Hz]})$'), 'interpreter', 'latex');
    ylabel(sprintf('$FFR(^{ENV}/_{TFS})$'), 'interpreter', 'latex');
    
    set(gca, 'fontsize', fSize);
    text(.1,.9,sprintf('p=%.4f, adj-R^2=%.4f', mdl.Coefficients.pValue(2), mdl.Rsquared.Adjusted), 'units', 'normalized', 'fontsize', fSize);
    grid on;
    
    set(gcf, 'units', 'normalized', 'position', [0.1 0.1 .8 .8]);
    if saveAgain
        saveas(gcf, [fig_save_dir figName], 'tiff');
    end
end

function curFilt= get_filter(fs_data)
N_bp_half= 4;
HalfPowerFrequency1=20;
HalfPowerFrequency2=2e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end