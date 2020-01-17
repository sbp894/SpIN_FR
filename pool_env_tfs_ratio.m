% Analyzing PTS vs NH pooled data of FFR(ENV/TFS) versus AUD(HF/LF)
clear;
clc;


CodesDirs= {'/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/helper', ...
    '/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/continuous'};
addpath(CodesDirs{:});


RootDataDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Artifact_Removed_FFR/';
CalibRootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
LatexDir= '/home/parida/Dropbox/Articles/Loss_of_tonotopy_in_HI_FFR/figures/';

saveData= 1;
saveFigs= 1;
mdl_scale_log0_lin1= 0;

fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/moving_segment_analysis_vowel/');
if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end

data_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Summary_slopes/';
if ~isfolder(data_save_dir)
    mkdir(data_save_dir);
end

allDirs= dir([RootDataDir '*SFR*']);
exlude_dirs= cell2mat(cellfun(@(x) contains(x,{'pink', 'q369_sfr_pilot1'}), lower({allDirs.name}'), 'uniformoutput', false));
allDirs= allDirs(~exlude_dirs);


subGroups.names= {'NH', 'PTS'};
subGroups.marker= {'ob', 'dr'};
subGroups.cols= {'b', 'r'};
subGroups.nums= length(subGroups.names);

band_cutoff_point= 500;

audio_freq_band_low= [60 band_cutoff_point];
audio_freq_band_high= [band_cutoff_point 5000];

data_f0_related_band= [60 band_cutoff_point];
warning('Assuming using artifact removed data');

t_latency= 5e-3;
useCalib= 0;
if useCalib
    baseSPL= 20*log10(1/sqrt(2)/(20e-6)); % For 1V wav-file, this is assumed dB SPL. Compare with calib file for correction.
    calibFigHan= 417;
    figure(calibFigHan); clf;
end

[sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');
restricted_time= find_voicing_boundaries(sig, fs, 0); %[0 .16; .2 .72; 0.88 1.2]; % exclude 50 ms onset

if ~isempty(restricted_time)
    warning('Using restricted time');
end

t_sig= (1:length(sig))/fs;
stim_dur= length(sig)/fs;

windowLength= 64e-3;
fracOverLap= 0;
nfft= 2^nextpow2(round(fs*windowLength));
NW= 1.5; % NW= dur * f_res => f_res= NW/dur. If NW=1.5, dur=50ms, f_res= 30 Hz
[~, PSDfreq] = pmtm(randn(round(fs*windowLength),1), NW, nfft, fs);
freq_inds2use= PSDfreq~=0;
PSDfreq=PSDfreq(freq_inds2use);


fracSlide= 1-fracOverLap;
tSlide= fracSlide*windowLength;
nSegs= 1 + floor((stim_dur-windowLength)/(tSlide));

pool_ratio_lf_to_hf_audio= nan(length(allDirs), nSegs);
pool_ratio_tfs_to_env_ffr= nan(length(allDirs), nSegs);
pool_hf_power_audio= nan(length(allDirs), nSegs);
pool_lf_power_audio= nan(length(allDirs), nSegs);
pool_env_power_ffr= nan(length(allDirs), nSegs);
pool_tfs_power_ffr= nan(length(allDirs), nSegs);

parfor dirVar= 1:length(allDirs)
% for dirVar= length(allDirs)
    curDir=  [RootDataDir allDirs(dirVar).name filesep];
    s_files= dir([curDir 'a*_S_*1*']); % choose dirs with attentuation=10 dB
    
    ratio_lf_to_hf_audio    = nan(nSegs, 1);
    ratio_tfs_to_env_f0     = nan(nSegs, 1);
    power_sig_low           = nan(nSegs, 1);
    power_sig_high          = nan(nSegs, 1);
    Pxx_env_seg             = nan(nSegs, 1);
    Pxx_tfs_seg             = nan(nSegs, 1);
    
    if ~isempty(s_files)
        
        if useCalib
            curCalibDir= [CalibRootDataDir allDirs(dirVar).name filesep];
            CalibFiles= dir([curCalibDir '*calib*']);
            if length(CalibFiles)~=1
                warning('Check if the correct calib file is choosen');
            end
            CalibFiles= CalibFiles(end).name;
            fprintf('Using %s as calib for file --> %s\n', CalibFiles, allDirs(dirVar).name);
            calibData= load([curCalibDir CalibFiles]);
            calibFreq= [0; calibData.data.CalibData(:, 1)*1e3; fs/2];
            calibGain= calibData.data.CalibData([1 1:end end], 2) - baseSPL;
            PSDgain= interp1(calibFreq, calibGain, PSDfreq);
            
            figure(calibFigHan);
            plot(calibData.data.CalibData(:, 1), calibData.data.CalibData(:, 2), 'DisplayName', allDirs(dirVar).name, 'LineWidth', 2);
            hold on;
        else 
            PSDgain= 0;
            fprintf('Using no calib for file --> %s\n', allDirs(dirVar).name);
        end
        
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
            if ~isempty(restricted_time)
                seg_in_valid_time= (seg_t_mid>restricted_time(:,1)) & (seg_t_mid<restricted_time(:,2));
            else
                seg_in_valid_time= 1;
            end
            
            if any(seg_in_valid_time)
                
                seg_stim= sig(seg_ind_start:seg_ind_end);
                seg_t= t_sig(seg_ind_start:seg_ind_end);
                
                cur_data_inds= t_data>(seg_t_start+t_latency) & t_data<(seg_t_end+t_latency);
                cur_data_env= s_data_env(cur_data_inds);
                cur_data_tfs= s_data_tfs(cur_data_inds);
                cur_t_data= t_data(cur_data_inds);
                
                
                [Pxx_sig_dB, freq_stim, ~]= plot_dpss_psd(seg_stim, fs, 'NW', NW, 'plot', false, 'nfft', nfft);
                Pxx_sig= 10.^((Pxx_sig_dB + PSDgain) / 10)*fs/nfft; % Can also use db2pow
                
                sig_freq_inds_low= freq_stim>audio_freq_band_low(1) & freq_stim<audio_freq_band_low(2);
                sig_freq_inds_high= freq_stim>audio_freq_band_high(1) & freq_stim<audio_freq_band_high(2);
                Pxx_sig_low= sum(Pxx_sig(sig_freq_inds_low)); 
                Pxx_sig_high= sum(Pxx_sig(sig_freq_inds_high));
                
                power_sig_low(segVar)= Pxx_sig_low;
                power_sig_high(segVar)= Pxx_sig_high;
                ratio_lf_to_hf_audio(segVar)= power_sig_low(segVar)/power_sig_high(segVar);
                
                
                Pxx_env_dB= plot_dpss_psd(cur_data_env, fs_data, 'NW', NW, 'plot', false); %#ok<*ASGLU>
                Pxx_env= 10.^(Pxx_env_dB/10); % Can also use db2pow
                
                [Pxx_tfs_dB, freq_data]= plot_dpss_psd(cur_data_tfs, fs_data, 'NW', NW, 'plot', false);
                Pxx_tfs= 10.^(Pxx_tfs_dB/10); % Can also use db2pow
                
                data_freq_inds= freq_data>data_f0_related_band(1) & freq_data<data_f0_related_band(2);
                Pxx_env_seg(segVar)= sum(Pxx_env(data_freq_inds));
                Pxx_tfs_seg(segVar)= sum(Pxx_tfs(data_freq_inds));
                
                ratio_tfs_to_env_f0(segVar)= Pxx_tfs_seg(segVar) / Pxx_env_seg(segVar);
            else
                
            end
        end
    end
    pool_ratio_lf_to_hf_audio(dirVar,:) = ratio_lf_to_hf_audio;
    pool_ratio_tfs_to_env_ffr(dirVar, :) = ratio_tfs_to_env_f0;
    pool_hf_power_audio(dirVar, :)      = power_sig_high;
    pool_lf_power_audio(dirVar, :)      = power_sig_low;
    pool_env_power_ffr(dirVar, :)       = Pxx_env_seg;
    pool_tfs_power_ffr(dirVar, :)       = Pxx_tfs_seg;
end

allChinData= repmat(struct('chinID', []), length(allDirs), 1);
for dirVar=1:length(allDirs)
    allChinData(dirVar).dirName= allDirs(dirVar).name;
    allChinData(dirVar).chinID= cell2mat(cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp(allDirs(dirVar).name,'(-Q\d+_)','tokens'), 'UniformOutput', 0));
    allChinData(dirVar).hf_power_audio= pool_hf_power_audio(dirVar, :);
    allChinData(dirVar).lf_power_audio= pool_lf_power_audio(dirVar, :);
    allChinData(dirVar).env_power_ffr= pool_env_power_ffr(dirVar, :);
    allChinData(dirVar).tfs_power_ffr= pool_tfs_power_ffr(dirVar, :);
    
    if contains(allChinData(dirVar).dirName, 'NH')
        allChinData(dirVar).group= 'NH';
    elseif contains(allChinData(dirVar).dirName, {'HI', 'PTS'})
        allChinData(dirVar).group= 'PTS';
    else
        warning('No group for %s', allChinData(dirVar).dirName);
        allChinData(dirVar).group= 'unknown';
    end
end


if saveData
    save([data_save_dir 'all_chins_data.mat'], 'allChinData');
end

%% plot
remove_outlier= true;
if remove_outlier
    warning('Debugging: excluding one PTS data');
    Exclude_point.chinID= [369];
    Exclude_point.type= 'PTS';
    Exclude_point.index= contains({allDirs.name}, Exclude_point.type) & contains({allDirs.name}, num2str(Exclude_point.chinID));
    pool_ratio_lf_to_hf_audio(Exclude_point.index, :)= nan;
    pool_ratio_tfs_to_env_ffr(Exclude_point.index, :)= nan;
    figure(2);
    figName= 'pooled_aud_env_ratios_outlier_excluded';
    calibStr= {'woCal', 'wCal'};
    figName_indv= sprintf('pooled_aud_env_ratios_%s', calibStr{useCalib+1});
else
    figure(1);
    figName= 'pooled_aud_env_ratios_all_included';
    figName_indv= sprintf('pooled_aud_env_ratios_%s', calibStr{useCalib+1});
end

clf;
co= get(gca, 'colororder');
markSize= 12;
fSize= 20;
tick_len= [.025 .025];
ax= nan(subGroups.nums, 1);
lw= 2;
lw2= 3;
xtick_val= [.001 .01 .1 1 10 100];
ytick_val= [.01 .1 1 10];
spLetters= 'ABCD';

matGreen= [0.4660, 0.6740, 0.1880];
matPurple= [0.4940, 0.1840, 0.5560];

% remove_outlier= true;
if mdl_scale_log0_lin1
    % leave as is
    
else
    pool_ratio_lf_to_hf_audio   = db(pool_ratio_lf_to_hf_audio);
    pool_ratio_tfs_to_env_ffr    = db(pool_ratio_tfs_to_env_ffr);
    xtick_val= db(xtick_val);
    ytick_val= db(ytick_val);
end

for typeVar= 1:subGroups.nums
    cur_type_inds= contains({allDirs.name}', subGroups.names(typeVar));
    ax(typeVar)= subplot( 1 , subGroups.nums, typeVar);
    cur_subgroup_data_x= pool_ratio_lf_to_hf_audio(cur_type_inds,:);
    cur_subgroup_data_y= pool_ratio_tfs_to_env_ffr(cur_type_inds,:);
    nanInds= isnan(cur_subgroup_data_x);
    
    cur_subgroup_data_x= cur_subgroup_data_x(~nanInds);
    cur_subgroup_data_y= cur_subgroup_data_y(~nanInds);
    
    cur_subgroup_est_x= sort(cur_subgroup_data_x);
    mdl= fitlm(cur_subgroup_data_x, cur_subgroup_data_y);
    c_m = mdl.Coefficients.Estimate;
    cur_subgroup_est_y= c_m(1)+ c_m(2)*cur_subgroup_est_x;
    
    plot(cur_subgroup_data_x, cur_subgroup_data_y, subGroups.marker{typeVar}, 'markersize', markSize, 'linew', lw);
    hold on;
    plot(cur_subgroup_est_x, cur_subgroup_est_y, '-', 'color', subGroups.cols{typeVar}, 'linew', lw2);
    
    grid off;
    
    if mdl_scale_log0_lin1
        set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fSize, 'box', 'off');
    else
        set(gca, 'fontsize', fSize);
    end
    title(strrep(subGroups.names{typeVar}, 'PTS', 'HI'));
    %     xlabel(sprintf('$Carrier( ^{HF_{power}[.5-3kHz]}/_{LF_{power}[50-500Hz]})$'), 'interpreter', 'latex');
    %     ylabel(sprintf('$FFR(^{ENV_{power}}/_{TFS_{power}})$'), 'interpreter', 'latex');
    %
    xlabel(sprintf('$LF_{stimulus}^{power}/HF_{stimulus}^{power}$'), 'interpreter', 'latex');
    if strcmp(subGroups.names(typeVar), 'NH')
        ylabel(sprintf('$TFS_{FFR}/ENV_{FFR}$'), 'interpreter', 'latex');
    end
    
    pValThresh= 1e-4;
    if mdl.Coefficients.pValue(2)>pValThresh
        text(.05,.9,sprintf('$p=%.3f, R^2=%.4f$', mdl.Coefficients.pValue(2), mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex');
    else
        text(.05,.9,sprintf('$p<%.3f, R^2=%.4f$', pValThresh, mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize, 'interpreter', 'latex');
    end
end

if remove_outlier
    %     title([subGroups.names{typeVar} '*']);
    title([strrep(subGroups.names{typeVar}, 'PTS', 'HI') '*']);
end

linkaxes(ax);
xlim([min(pool_ratio_lf_to_hf_audio(:))*.9 max(pool_ratio_lf_to_hf_audio(:))*1.1]);
ylim([min(pool_ratio_tfs_to_env_ffr(:))*.9 max(pool_ratio_tfs_to_env_ffr(:))*1.1]);

set(gcf, 'units', 'normalized', 'position', [0.1 0.1 .8 .8]);

saveas(gcf, [fig_save_dir figName], 'png');
if saveFigs
    saveas(gcf, [LatexDir figName], 'epsc');
end

%% plot individually
figure(6156);
clf;
lw3= 5;
xtick_label= cellfun(@(x) num2str(x), num2cell(xtick_val), 'uniformoutput', false);
ytick_label= cellfun(@(x) num2str(x), num2cell(ytick_val), 'UniformOutput', false);
fSize_txt= 20;
fSize= 20;

subplot( 1 , subGroups.nums, 1);
co= set_colblind_order();
subplot( 1 , subGroups.nums, 2);
set_colblind_order();

for typeVar= 1:subGroups.nums
    cur_type_inds= find(contains({allDirs.name}', subGroups.names(typeVar)));
    
    if strcmp(subGroups.names{typeVar}, 'NH')
        col_val= co(1,:);
    else
        col_val= co(2,:);
%         col_val= get_colormap([1 0.1 0.1], [0.8 0.1 0.1], length(cur_type_inds));
    end
    
    %     set(gca, 'colororder', col_ord);
    
    ax(typeVar)= subplot( 1 , subGroups.nums, typeVar);
    
    for chinVar = 1:length(cur_type_inds)
        cur_subgroup_data_x= pool_ratio_lf_to_hf_audio(cur_type_inds(chinVar),:);
        cur_subgroup_data_y= pool_ratio_tfs_to_env_ffr(cur_type_inds(chinVar),:);
        
        nanInds= isnan(cur_subgroup_data_x);
        
        cur_subgroup_data_x= cur_subgroup_data_x(~nanInds);
        cur_subgroup_data_y= cur_subgroup_data_y(~nanInds);
        if ~isempty(cur_subgroup_data_x)
            
            cur_subgroup_est_x= sort(cur_subgroup_data_x);
            mdl= fitlm(cur_subgroup_data_x, cur_subgroup_data_y);
            c_m = mdl.Coefficients.Estimate;
            cur_subgroup_est_y= c_m(1)+ c_m(2)*cur_subgroup_est_x;
            
            plot(cur_subgroup_data_x, cur_subgroup_data_y, subGroups.marker{typeVar}, 'markersize', markSize, 'linew', lw, 'color', col_val);
            hold on;
            plot(cur_subgroup_est_x, cur_subgroup_est_y, '-', 'color', col_val, 'linew', lw2);
            %             set(gca, 'XColor', matGreen, 'YColor', matPurple);
        end
    end
    
    grid off;
    if mdl_scale_log0_lin1
        set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', fSize, 'XTick', xtick_val, 'XTickLabel', xtick_label, ...
            'YTick', ytick_val, 'YTickLabel', ytick_label, 'box', 'off', 'LineWidth', 1.5, 'TickLength', tick_len);
    else
        set(gca, 'fontsize', fSize, 'XTick', xtick_val, 'XTickLabel', xtick_label, 'YTick', ytick_val, ...
            'YTickLabel', ytick_label, 'box', 'off', 'LineWidth', 1.5, 'TickLength', tick_len);
    end
    title(strrep(subGroups.names{typeVar}, 'PTS', 'HI'));
    %     xlabel(sprintf('$Carrier( ^{HF_{power}[.5-3kHz]}/_{LF_{power}[50-500Hz]})$'), 'interpreter', 'latex');
    %     xlabel(sprintf('$Carrier( ^{HF_{power}}/_{LF_{power}})$'), 'interpreter', 'latex');
    %     ylabel(sprintf('$FFR(^{ENV_{power}}/_{TFS_{power}})$'), 'interpreter', 'latex');
    xlabel(sprintf('$LF_{stimulus}^{power}/HF_{stimulus}^{power}$ (dB)'), 'interpreter', 'latex');
    if strcmp(subGroups.names(typeVar), 'NH')
        ylabel(sprintf('$TFS_{FFR}^{power}/ENV_{FFR}^{power}$ (dB)'), 'interpreter', 'latex');
    end
    
    
    cur_subgroup_data_x= pool_ratio_lf_to_hf_audio(cur_type_inds,:);
    cur_subgroup_data_y= pool_ratio_tfs_to_env_ffr(cur_type_inds,:);
    nanInds= isnan(cur_subgroup_data_x);
    cur_subgroup_data_x= cur_subgroup_data_x(~nanInds);
    cur_subgroup_data_y= cur_subgroup_data_y(~nanInds);
    cur_subgroup_est_x= sort(cur_subgroup_data_x);
    mdl= fitlm(cur_subgroup_data_x, cur_subgroup_data_y);
    c_m = mdl.Coefficients.Estimate;
    cur_subgroup_est_y= c_m(1)+ c_m(2)*cur_subgroup_est_x;
    %     plot(cur_subgroup_est_x, cur_subgroup_est_y, '-', 'color', subGroups.cols{typeVar}, 'linew', lw3);
    plot(cur_subgroup_est_x, cur_subgroup_est_y, '-', 'color', 'k', 'linew', lw3);
    pValThresh= 1e-4;
    txtGap= .07;
    if mdl.Coefficients.pValue(2)>pValThresh
%         text(.05,.95,sprintf('$r=%.2f$', mdl.Coefficients.Estimate(2)), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold');
        text(.05,.95-txtGap,sprintf('$p=%.2f$', mdl.Coefficients.pValue(2)), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold');
        text(.05,.95-2*txtGap,sprintf('$R^2=%.2f$', mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold');
    else
%         text(.05,.95,sprintf('$r=%.2f$', mdl.Coefficients.Estimate(2)), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold');
        text(.05,.95-txtGap,sprintf('$p<%.4f$', pValThresh), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold');
        text(.05,.95-2*txtGap,sprintf('$R^2=%.2f$', mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold');
    end
end

strrep(subGroups.names{typeVar}, 'PTS', 'HI')

linkaxes(ax);
set(gcf, 'units', 'inches', 'position', [1 1 11 5]);
txtHan= add_subplot_letter(1, 2, 30);

saveas(gcf, [fig_save_dir figName_indv], 'png');
if saveFigs
    saveas(gcf, [LatexDir figName_indv], 'epsc');
end

nhInds= find(contains({allDirs.name}', 'NH'));
hiInds= find(contains({allDirs.name}', 'PTS') & ~Exclude_point.index');

nh_ChinIDs=cell2mat(cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp({allDirs(nhInds).name}','(-Q\d+_)','tokens'), 'UniformOutput', 0));
hi_ChinIDs=cell2mat(cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp({allDirs(hiInds).name}','(-Q\d+_)','tokens'), 'UniformOutput', 0));

test_env_tfs_audPower_relation...
    (pool_lf_power_audio, pool_hf_power_audio, pool_env_power_ffr, pool_tfs_power_ffr, nhInds, hiInds, LatexDir, calibStr{useCalib+1}, data_save_dir, nh_ChinIDs, hi_ChinIDs);
% test_env_tfs_audPower_relation_trimmer(pool_lf_power_audio, pool_hf_power_audio, pool_env_power_ffr, pool_tfs_power_ffr, nhInds, hiInds, LatexDir);
check_variability_in_mdls(pool_ratio_lf_to_hf_audio, pool_ratio_tfs_to_env_ffr, pool_lf_power_audio, pool_hf_power_audio, pool_env_power_ffr, pool_tfs_power_ffr, nhInds, hiInds)



rmpath(CodesDirs{:});


%%
function curFilt= get_filter(fs_data)
N_bp_half= 2;
HalfPowerFrequency1=80;
HalfPowerFrequency2=1e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end