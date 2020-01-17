clear;
clc;

% AllChinIDs= [358 360 365 366 367 368 369 370];
% AllChinIDs= [369 370]; % [366 369];
% AllChinIDs= 370;
CodesDirs= {'/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/helper', ...
    '/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/continuous'};
addpath(CodesDirs{:});
forceSameSegment_NH_HI= 0;


chin_groups= {'nh','pts'};
for groupVar= 1:length(chin_groups)
    if strcmpi(chin_groups{groupVar}, 'NH')
        % NH
        %         AllChinIDs= [365 368 369 370 374];
        AllChinIDs= 370;
    elseif strcmpi(chin_groups{groupVar}, 'PTS')
        % PTS
        AllChinIDs= 367; %[358 360 366 367 369 370];
    end
    
    mdl_scale_log0_lin1= 0;
    nMax2Plot= 1;
    
    [sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');
    t_sig= (1:length(sig))/fs;
    stim_dur= length(sig)/fs;
    % LatexDir= '/home/parida/Dropbox/Seminars/SHRP_Feb19/figures/';
%     LatexDir= '/home/parida/Dropbox/Conferences/ASA-2019/slides/figures/';
    LatexDir= '/home/parida/Dropbox/Articles/Loss_of_tonotopy_in_HI_FFR/figures/';
    
%     restricted_time= [0 .72; 0.88 1.2];
    restricted_time= find_voicing_boundaries(sig, fs, 0);
    % restricted_time= [0 1.3];
    saveAgain= 1;
    
    windowLength= 64e-3;
    fracOverLap= 0;
    NW= 1.5; % NW= dur * f_res => f_res= NW/dur. If NW=1.5, dur=50ms, f_res= 30 Hz
    
    fracSlide= 1-fracOverLap;
    tSlide= fracSlide*windowLength;
    nSegs= 1 + floor((stim_dur-windowLength)/(tSlide));
    
    dpss_yRange= 60;
    plot_dpss_each_iter= false;
    lw= 2;
    lw2= 3;
    nSProws=2;
    fSize= 20;
    tick_len1= [.015 .015];
    tick_len2= [.025 .025];
    nSPcols=1;
    matGreen= [0.4660, 0.6740, 0.1880];
    matPurple= [0.4940, 0.1840, 0.5560];
    
    t_latency= 5e-3;
    
    ratio_lf_to_hf_audio= nan(nSegs, 1);
    audio_freq_band_low= [60 500];
    audio_freq_band_high= [500 5000];
    
    data_f0_related_band= [60 500];
    ratio_f0_tfs_to_env= nan(nSegs, 1);
    power_f0_env= nan(nSegs, length(AllChinIDs));
    power_f0_tfs= nan(nSegs, length(AllChinIDs));
    
    slope_vals= nan(length(AllChinIDs), 1);
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
        end
        allFiles= allFiles( ~contains(lower({allFiles.name}'), 'pink') | ismember(chinID, [373 374])); % data from Q373 374 named as *pink*, but used combine_sfr_pics to create a file to use here
        if length(allFiles)>1
            
            fprintf('there are multiple dirs. \n');
            
            for dirVar= 1:length(allFiles)
                fprintf('(%d)-%s\n', dirVar, allFiles(dirVar).name);
            end
            if strcmpi(chin_groups{groupVar}, 'NH')
                chosen_dir_num=length(allFiles)-1;
            elseif strcmpi(chin_groups{groupVar}, 'PTS')
                chosen_dir_num=length(allFiles);
            end
            %             chosen_dir_num= input('Which one? \n');
            if contains(allFiles(chosen_dir_num).name, 'NH')
                postFix= 'NH';
            elseif contains(allFiles(chosen_dir_num).name, 'PTS')
                postFix= 'PTS';
            else
                postFix= 'NH';
                warning('Assuming NH for %s', allFiles(chosen_dir_num).name);
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
        fprintf('---------- Working on %s ----------\n', data_dir);
        
        fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/moving_segment_analysis_vowel/');
        if ~isfolder(fig_save_dir)
            mkdir(fig_save_dir);
        end
        art_fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/Artifact_removed/');
        if ~isfolder(art_fig_save_dir)
            mkdir(art_fig_save_dir);
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
        
        gain= 20e3/10; % divide by 10 because we thought dagan headstage has 10 gain, which is not included in the dagan knob. But actually it is.
        s_data_pos= s_data_pos(:)/gain*1e6; % now in microvolt
        s_data_neg= s_data_neg(:)/gain*1e6;
        
        initialRampDur= 20e-3;
        ramp_nSamples= round(initialRampDur*fs_data);
        rampHamming= hamming(2*ramp_nSamples)';
        rampVector= [rampHamming(1:ramp_nSamples), ones(1, length(s_data_pos)-length(rampHamming)) rampHamming(ramp_nSamples+1:end)]';
        curFilt= get_filter(fs_data);
        
        
        % --------------------------
        % Remove 75 Hz and its harmonics
        plotVar= 1;
        if plotVar
            %             warning('debugging');
        end
        artifact_fig= 44;
        figure(artifact_fig); clf;
        subplot(211);
        s_data_pos=remove_artifact_ffr(s_data_pos, fs_data, plotVar);
        set(gca, 'FontSize', fSize);
        title('$FFR_{+ve}$', 'Interpreter', 'latex');
        subplot(212);
        s_data_neg=remove_artifact_ffr(s_data_neg, fs_data, plotVar);
        set(gca, 'FontSize', fSize);
        title('$FFR_{-ve}$', 'Interpreter', 'latex');
        legend('Raw', 'Processed', 'Location', 'southwest');
        set(artifact_fig, 'Units', 'inches', 'Position', [1 1 9 8]);
        
        ArtFigName= sprintf('Q%d_atten%.0fdB_%s_ArtRemove', chinID, s_atten, postFix);
        saveas(gcf, [art_fig_save_dir ArtFigName], 'png');
        close(artifact_fig);
        % --------------------------
        s_data_pos_filt= filtfilt(curFilt, s_data_pos.*rampVector).*rampVector;
        s_data_neg_filt= filtfilt(curFilt, s_data_neg.*rampVector).*rampVector;
        s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
        s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;
        t_data= (1:length(s_data_env))/fs_data;
        
        A_stim= max([s_data_env;s_data_tfs]);
        A_shift= 1.3*A_stim;
        
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
                ratio_lf_to_hf_audio(segVar)= Pxx_sig_low/Pxx_sig_high;
                
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
                
                ratio_f0_tfs_to_env(segVar)= Pxx_tfs_seg/Pxx_env_seg;
                power_f0_env(segVar, chinVar)= Pxx_env_seg;
                power_f0_tfs(segVar, chinVar)= Pxx_tfs_seg;
                
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
        %%y
        xtick_vals_ratio= [.1 1 1e1 1e2 1e3];
        ytick_vals_ratio= [.01 .1 1 10];
        if mdl_scale_log0_lin1
            % leave as is
            
        else
            ratio_lf_to_hf_audio= db(ratio_lf_to_hf_audio);
            ratio_f0_tfs_to_env= db(ratio_f0_tfs_to_env);
            xtick_vals_ratio= db(xtick_vals_ratio);
            ytick_vals_ratio= db(ytick_vals_ratio);
        end
        
        
        ratio_lf_to_hf_audio_est= sort(ratio_lf_to_hf_audio);
        mdl= fitlm(ratio_lf_to_hf_audio, ratio_f0_tfs_to_env);
        c_m = mdl.Coefficients.Estimate;
        ratio_f0_tfs_to_env_est= c_m(1)+ c_m(2)*ratio_lf_to_hf_audio_est;
        slope_vals(chinVar)= c_m(2);
        
        if ~forceSameSegment_NH_HI
            [~, ratio_f0_tfs_to_env_max_inds]= sort(ratio_f0_tfs_to_env, 'descend','MissingPlacement','last');
%             [~, ratio_f0_tfs_to_env_max_inds]= sort(ratio_f0_tfs_to_env, 'ascend','MissingPlacement','last');
            ratio_f0_tfs_to_env_max_inds= ratio_f0_tfs_to_env_max_inds(1:nMax2Plot); % Plot only nMax2Plot segs
        else
            ratio_f0_tfs_to_env_max_inds = 17;
        end
        seg_ind_start= nan(nMax2Plot, 1);
        seg_ind_end= nan(nMax2Plot, 1);
        for maxVar=1:length(ratio_f0_tfs_to_env_max_inds)
            seg_t_start=  (ratio_f0_tfs_to_env_max_inds(maxVar)-1)*tSlide;
            seg_ind_start(maxVar)= max(1, round(seg_t_start*fs));
            seg_t_end= seg_t_start + windowLength;
            seg_ind_end(maxVar)= min(length(sig), round(seg_t_end*fs));
        end
        
        figure(2);
        clf;
        subplot(211);
        
        hold on
        plot(t_sig, -A_shift + A_stim*sig, 'color', 'k'); % very first time
        for maxVar=1:nMax2Plot
            plot(t_sig(seg_ind_start(maxVar):seg_ind_end(maxVar)), -A_shift + A_stim*sig(seg_ind_start(maxVar):seg_ind_end(maxVar)), 'r'); % very first time
            plot( [t_sig(seg_ind_start(maxVar)) t_sig(seg_ind_start(maxVar))], [-2*A_shift  2*A_shift], 'r', 'linew', lw);
            plot( [t_sig(seg_ind_end(maxVar)) t_sig(seg_ind_end(maxVar))], [-2*A_shift  2*A_shift], 'r', 'linew', lw);
%             text(t_sig(round((seg_ind_start(maxVar)+seg_ind_end(maxVar))/2)), -1.1*A_shift + A_stim*sig(round((seg_ind_start(maxVar)+seg_ind_end(maxVar))/2)), sprintf('%d', maxVar));
        end
        plot(t_data, A_shift + s_data_env, 'color', get_color('purple')); % very first time
        plot(t_data, s_data_tfs, 'color', get_color('dg')); % very first time
%         title(sprintf('Animal-ID:%d, SPL=%.0f dB, %s', chinID, 75-s_atten, postFix));
        title('');
        
        lg_hans(3)= plot(nan, nan, 'color', 'k', 'linew', lw);
        lg_hans(1)= plot(nan, nan, 'color', get_color('purple'), 'linew', lw);
        lg_hans(2)= plot(nan, nan, 'color', get_color('dg'), 'linew', lw);
        lg= legend(lg_hans, 'ENV', 'TFS', 'SIG');
        grid off;
        
        figName= sprintf('Q%d_atten%.0fdB_%s_nSeg%d', chinID, s_atten, postFix, nMax2Plot);
        set(gca, 'fontsize', fSize, 'LineWidth', 1.5, 'TickLength', tick_len1);
        xlabel('Time (sec)');
        set(lg, 'fontsize', fSize, 'location', 'southeast', 'box', 'off');
        axis tight;
        ylabel('FFR-Amp (\muV)');
        text(0, 1.075, '(\bfA\rm)', 'FontSize', 30, 'Units', 'normalized');
        
        subplot(223);
        %     yyaxis left;
        hold on;
        fill([.5 11 11 .5 .5]*1e3, [-100 -100 0 0 -150], [.8 .8 .9], 'FaceAlpha', .25);
        [Pxx_sig, freq_stim, px_stim]= plot_dpss_psd(sig, fs, 'NW', NW);
        ylim([max(Pxx_sig)+10-dpss_yRange max(Pxx_sig)+10])
        ylabel('Stim-PSD (dB/Hz)');
        set(px_stim, 'Color', 'r');
        
        %     yyaxis right;
        %     hold on;
        %     freq_ticks= [1 10 100 1e3];
        %     [Pxx_env, ~, px_env]= plot_dpss_psd(s_data_env, fs_data, 'NW', NW);
        %     set(px_env, 'color', co(4,:), 'linestyle', '-');
        %     [Pxx_tfs, freq_data, px_tfs]= plot_dpss_psd(s_data_tfs, fs_data, 'NW', NW);
        %     ylim([max([Pxx_env;Pxx_tfs])+10-dpss_yRange max([Pxx_env;Pxx_tfs])+10])
        %     set(px_tfs, 'color', co(5,:), 'linestyle', '-');
        %     lg= legend('SIG', 'ENV', 'TFS');
        %     set(lg, 'fontsize', 12, 'location', 'northwest');
        %     set(gca, 'fontsize', fSize, 'ytick', [], 'xtick', freq_ticks);
        %     ylabel('');
        %     title('');
        xlabel('Frequency (Hz)');
        xtick_vals= [50 500 5e3];
        xtick_labs= cellfun(@(x) num2str(x), num2cell(xtick_vals), 'UniformOutput', false);
        set(gca, 'XTick', xtick_vals, 'XTickLabel', xtick_labs, 'LineWidth', 1.5, 'TickLength', tick_len2);
        xlim([.99*min(xtick_vals) 1.01*max(xtick_vals)]);
        set(gca, 'fontsize', fSize);
        title('');
        grid off;
        
        subplot(224);
        mrkSize= 16;
        hold on;
        plot(ratio_lf_to_hf_audio, ratio_f0_tfs_to_env, 'o', 'color', get_color('gray'), 'markersize', mrkSize, 'linew', lw);
        for maxVar=1:nMax2Plot
            plot(ratio_lf_to_hf_audio(ratio_f0_tfs_to_env_max_inds(maxVar)), ratio_f0_tfs_to_env(ratio_f0_tfs_to_env_max_inds(maxVar)), 'ro', 'markersize', mrkSize, 'linew', lw2);
%             text(ratio_lf_to_hf_audio(ratio_f0_tfs_to_env_max_inds(maxVar)), ratio_f0_tfs_to_env(ratio_f0_tfs_to_env_max_inds(maxVar)), sprintf('%d', maxVar));
        end
        plot(ratio_lf_to_hf_audio_est, ratio_f0_tfs_to_env_est, '-k', 'linew', lw2);
        ytick_labs= cellfun(@(x) num2str(x), num2cell(ytick_vals_ratio), 'UniformOutput', false);
        
        if mdl_scale_log0_lin1
            set(gca, 'xscale', 'log', 'yscale', 'log', 'ytick', ytick_vals_ratio, 'YTickLabel',ytick_labs, 'LineWidth', 1.5, 'TickLength', tick_len2);
        else
            set(gca, 'ytick', ytick_vals_ratio, 'YTickLabel',ytick_labs, 'LineWidth', 1.5, 'TickLength', tick_len2);
        end
        
        ylim([min(ytick_vals_ratio) max(ytick_vals_ratio)]);
        
        subplot(223);
        text(0, 1.075, '(\bfB\rm)', 'FontSize', 30, 'Units', 'normalized');
        subplot(224);
        text(0, 1.075, '(\bfC\rm)', 'FontSize', 30, 'Units', 'normalized');
        
        
%         set(gca, 'XColor', matGreen, 'YColor', matPurple);
        
        %     xlabel(sprintf('Audio ratio (norm) HF (%.1f-%.1f kHz) to LF (%.2f-%.2f kHz)', audio_freq_band_high(1)/1e3, audio_freq_band_high(2)/1e3, audio_freq_band_low(1)/1e3, audio_freq_band_low(2)/1e3));
        %     ylabel(sprintf('ENV/TFS in FFR in %.2f-%.2f kHz', data_f0_related_band(1)/1e3, data_f0_related_band(2)/1e3));
        
        %     xlabel(sprintf('$Stimulus( ^{HF_{power}[.5-3kHz]}/_{LF_{power}[50-500Hz]})$'), 'interpreter', 'latex');
        %     xlabel(sprintf('$Stimulus( ^{HF_{power}}/_{LF_{power}})$'), 'interpreter', 'latex');
        %     ylabel(sprintf('$FFR(^{ENV_{power}}/TFS_{power}$'), 'interpreter', 'latex');
        xlabel(sprintf('$LF_{stimulus}^{power}/HF_{stimulus}^{power}$ (dB)'), 'interpreter', 'latex');
        ylabel(sprintf('$TFS_{FFR}^{power}/ENV_{FFR}^{power}$ (dB)'), 'interpreter', 'latex');
        
        xtick_labs= cellfun(@(x) num2str(x), num2cell(xtick_vals_ratio), 'UniformOutput', false);
        set(gca, 'fontsize', fSize, 'xtick', xtick_vals_ratio, 'XTickLabel', xtick_labs);
        pValThresh= 1e-3;
        if mdl.Coefficients.pValue(2)>1e-3
%             text(.1,.1,sprintf('$r=%.2f, p=%.2f, R^2=%.2f$', mdl.Coefficients.Estimate(2), mdl.Coefficients.pValue(2), mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex');
            text(.1,.1,sprintf('$p=%.2f, R^2=%.2f$', mdl.Coefficients.pValue(2), mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex');
        else
%             text(.1,.1,sprintf('$r=%.2f, p<%.3f, R^2=%.2f$', mdl.Coefficients.Estimate(2), pValThresh, mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex');
            text(.1,.1,sprintf('$p<%.3f, R^2=%.2f$', pValThresh, mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex');
        end
        
        grid off
        %     xlim([1e-3 2]);
        
        set(gcf, 'units', 'normalized', 'position', [0.1 0.1 .8 .8]);
        
        
        if saveAgain
            saveas(gcf, [fig_save_dir figName], 'png');
            if chinID==370 && strcmp(postFix, 'NH') && nMax2Plot==1
                saveas(gcf, [LatexDir figName], 'epsc');
            end
            if chinID==367 && strcmp(postFix, 'PTS') && nMax2Plot==1
                saveas(gcf, [LatexDir figName], 'epsc');
            end
            if nMax2Plot==0
                saveas(gcf, [LatexDir figName], 'epsc');
            end
        end
    end
end
rmpath(CodesDirs{:});

function curFilt= get_filter(fs_data)
N_bp_half= 2;
HalfPowerFrequency1=70;
HalfPowerFrequency2=1e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end