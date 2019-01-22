% Analyzing PTS vs NH pooled data of FFR(ENV/TFS) versus AUD(HF/LF)
clear;
clc;

RootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';

fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/moving_segment_analysis_vowel/');
if ~isdir(fig_save_dir)
    mkdir(fig_save_dir);
end

allDirs= dir([RootDataDir '*SFR*']);
exlude_dirs= cell2mat(cellfun(@(x) contains(x,'pink'), lower({allDirs.name}'), 'uniformoutput', false));
allDirs= allDirs(~exlude_dirs);


subGroups.names= {'NH', 'PTS'};
subGroups.marker= {'ob', 'dr'};
subGroups.cols= {'b', 'r'};
subGroups.nums= length(subGroups.names);
% subGroups.data= repmat(struct('ratio_hf_to_lf_audio', [], 'ratio_f0_env_to_tfs', []), length(allDirs), 1);

audio_freq_band_low= [50 500];
audio_freq_band_high= [500 3000];

data_f0_related_band= [50 500];

restricted_time= [0 .72; 0.88 1.2];
t_latency= 5e-3;
% restricted_time= [0 1.3];

[sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');
t_sig= (1:length(sig))/fs;
stim_dur= length(sig)/fs;

windowLength= 65e-3;
fracOverLap= .75;
NW= 1.5; % NW= dur * f_res => f_res= NW/dur. If NW=1.5, dur=50ms, f_res= 30 Hz

fracSlide= 1-fracOverLap;
tSlide= fracSlide*windowLength;
nSegs= 1 + floor((stim_dur-windowLength)/(tSlide));

pool_ratio_hf_to_lf_audio= nan(length(allDirs), nSegs);
pool_ratio_f0_env_to_tfs= nan(length(allDirs), nSegs);

parfor dirVar= 1:length(allDirs)
    curDir=  [RootDataDir allDirs(dirVar).name filesep];
    s_files= dir([curDir 'a*_S_*10*']); % choose dirs with attentuation=10 dB
    
    ratio_hf_to_lf_audio= nan(nSegs, 1);
    ratio_f0_env_to_tfs= nan(nSegs, 1);
    
    if ~isempty(s_files)
        fprintf('------- Working for %s \n', allDirs(dirVar).name);
        
        s_data_cell= cell(length(s_files), 2);
        nPairs_actual= nan(length(s_files), 1);
        
        for sfile_var=1:length(s_files)
            temp_data= load([curDir s_files(sfile_var).name]);
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
                
                
                [Pxx_sig_dB, freq_stim, ~]= plot_dpss_psd(seg_stim, fs, 'NW', NW, 'plot', false);
                Pxx_sig= 10.^(Pxx_sig_dB/10);
                
                sig_freq_inds_low= freq_stim>audio_freq_band_low(1) & freq_stim<audio_freq_band_low(2);
                sig_freq_inds_high= freq_stim>audio_freq_band_high(1) & freq_stim<audio_freq_band_high(2);
                Pxx_sig_low= sum(Pxx_sig(sig_freq_inds_low));
                Pxx_sig_high= sum(Pxx_sig(sig_freq_inds_high));
                ratio_hf_to_lf_audio(segVar)= (Pxx_sig_high/range(audio_freq_band_high))/(Pxx_sig_low/range(audio_freq_band_low));
                
                
                Pxx_env_dB= plot_dpss_psd(cur_data_env, fs_data, 'NW', NW, 'plot', false); %#ok<*ASGLU>
                Pxx_env= 10.^(Pxx_env_dB/10);
                
                [Pxx_tfs_dB, freq_data]= plot_dpss_psd(cur_data_tfs, fs_data, 'NW', NW, 'plot', false);
                Pxx_tfs= 10.^(Pxx_tfs_dB/10);
                
                data_freq_inds= freq_data>data_f0_related_band(1) & freq_data<data_f0_related_band(2);
                Pxx_env_seg= sum(Pxx_env(data_freq_inds));
                Pxx_tfs_seg= sum(Pxx_tfs(data_freq_inds));
                
                ratio_f0_env_to_tfs(segVar)= Pxx_env_seg/Pxx_tfs_seg;
                
            else
                ratio_hf_to_lf_audio(segVar)= nan;
                ratio_f0_env_to_tfs(segVar)= nan;
            end
        end
    end
    pool_ratio_hf_to_lf_audio(dirVar, :)= ratio_hf_to_lf_audio;
    pool_ratio_f0_env_to_tfs(dirVar, :)= ratio_f0_env_to_tfs;
end


%% plot
remove_outlier= true;
if remove_outlier
    warning('Debugging: excluding one PTS data');
    Exclude_point.chinID= 369;
    Exclude_point.type= 'PTS';
    Exclude_point.index= contains({allDirs.name}, Exclude_point.type) & contains({allDirs.name}, num2str(Exclude_point.chinID));
    pool_ratio_hf_to_lf_audio(Exclude_point.index, :)= nan;
    pool_ratio_f0_env_to_tfs(Exclude_point.index, :)= nan;
    figure(2);
    figName= 'pooled_aud_env_ratios_all_included';
else
    figure(1);
    figName= 'pooled_aud_env_ratios_outlier_excluded';
end

clf;
markSize= 12;
fSize= 14;
ax= nan(subGroups.nums, 1);
lw= 3;

% remove_outlier= true;


for typeVar= 1:subGroups.nums
    cur_type_inds= contains({allDirs.name}', subGroups.names(typeVar));
    ax(typeVar)= subplot( 1 , subGroups.nums, typeVar);
    cur_subgroup_data_x= pool_ratio_hf_to_lf_audio(cur_type_inds,:);
    cur_subgroup_data_y= pool_ratio_f0_env_to_tfs(cur_type_inds,:);
    nanInds= isnan(cur_subgroup_data_x);
    
    cur_subgroup_data_x= cur_subgroup_data_x(~nanInds);
    cur_subgroup_data_y= cur_subgroup_data_y(~nanInds);
    
    cur_subgroup_est_x= sort(cur_subgroup_data_x);
    mdl= fitlm(cur_subgroup_data_x, cur_subgroup_data_y);
    c_m = mdl.Coefficients.Estimate;
    cur_subgroup_est_y= c_m(1)+ c_m(2)*cur_subgroup_est_x;
    
    plot(cur_subgroup_data_x, cur_subgroup_data_y, subGroups.marker{typeVar}, 'markersize', markSize);
    hold on;
    plot(cur_subgroup_est_x, cur_subgroup_est_y, '-', 'color', subGroups.cols{typeVar}, 'linew', lw);
    
    grid on;
    
    set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fSize);
    title(subGroups.names{typeVar});
    xlabel(sprintf('Audio ratio (norm) HF (%.1f-%.1f kHz) to LF (%.2f-%.2f kHz)', audio_freq_band_high(1)/1e3, audio_freq_band_high(2)/1e3, audio_freq_band_low(1)/1e3, audio_freq_band_low(2)/1e3));
    ylabel(sprintf('ENV/TFS in FFR in %.2f-%.2f kHz', data_f0_related_band(1)/1e3, data_f0_related_band(2)/1e3));
    text(.1,.9,sprintf('p=%.4f, adj-R^2=%.4f', mdl.Coefficients.pValue(2), mdl.Rsquared.Adjusted), 'units', 'normalized', 'fontsize', fSize);
    
end
if remove_outlier
    title([subGroups.names{typeVar} '(removed one outlier (Q369, PTS))']);
end

linkaxes(ax);
xlim([min(pool_ratio_hf_to_lf_audio(:))*.9 max(pool_ratio_hf_to_lf_audio(:))*1.1]);
ylim([min(pool_ratio_f0_env_to_tfs(:))*.9 max(pool_ratio_f0_env_to_tfs(:))*1.1]);

set(gcf, 'units', 'normalized', 'position', [0.1 0.1 .8 .8]);
saveas(gcf, [fig_save_dir figName], 'tiff');
%%
function curFilt= get_filter(fs_data)
N_bp_half= 4;
HalfPowerFrequency1=20;
HalfPowerFrequency2=2e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end