% tests effect of tdt-generated pink noise on SFR

clear;
clc;
figure(11);
clf;

allChinIDs= [371 373 374 379]; %[371 373 374];
RootDataDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Artifact_Removed_FFR/';
fig_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/tdt_pink_SFR/';
data_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/SFRpink500/';
LatexDir= '/home/parida/Dropbox/Articles/Loss_of_tonotopy_in_HI_FFR/figures/';
stimDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/stim/';

mrkOrder= 'ovdph';
tStart= 0; tEnd= 1.3; tNF= 1.3;
tStart_whole= 0; tEnd_whole= 1.3;
saveFigs= 1;
saveLatex= 1;
saveFinalFig= 1;
saveData= 0;
plotNF= 0;
remove_artifact_here= 0; % because using cleaned data => already artifact removed
if ~remove_artifact_here
    fprintf('Assuming using artifact removed Data\n');
end

if ~isfolder(data_save_dir)
    mkdir(data_save_dir);
end

if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end
if plotNF
    l1= nan(length(allChinIDs), 3);
    legStr= cell(length(allChinIDs), 3);
else
    l1= nan(length(allChinIDs), 2);
end
for chinVar= 1:length(allChinIDs)
    chinID= allChinIDs(chinVar);
    
    data_dir= dir([RootDataDir '*Q' num2str(chinID) '*pink500*']);
    data_dir= [RootDataDir data_dir.name filesep];
    
    allfiles= dir([data_dir 'a*SFR*.mat']);
    allfiles= allfiles(~(contains({allfiles.name}, 'latency') | contains({allfiles.name}, 'artifact')));
    
    all_snrs= cell2mat(cellfun(@(x) str2double(strrep(x(regexp(x, 'snr_')+4 : regexp(x, '_atn')-1), 'm', '-')), {allfiles.name}, 'uniformoutput', false));
    all_snrs(isnan(all_snrs))= [];
    all_snrs= fliplr(unique(all_snrs));
    
    fig_save_dir_subdir= [fig_save_dir sprintf('t%.0fto%.0f_ms/', tStart*1e3, tEnd*1e3)];
    
    
    raw_power_env= nan(length(all_snrs), 1);
    frac_power_env= nan(length(all_snrs), 1);
    raw_power_tfs= nan(length(all_snrs), 1);
    frac_power_tfs= nan(length(all_snrs), 1);
    % many way to compute NF. Get independent NF-PSD estimates for ENV and TFS
    % from different SNRs. Avereage them to obtain one NF for TFS and one for
    % ENV.
    raw_power_nf_env= nan(length(all_snrs), 1);
    raw_power_nf_tfs= nan(length(all_snrs), 1);
    
    raw_power_nf_env_CI= nan(length(all_snrs), 2);
    raw_power_nf_tfs_CI= nan(length(all_snrs), 2);
    
    stim_fName= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2017_11_02-Q325_AN_NH/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav';
    [sig, fs_sig]= audioread(stim_fName);
    
    snrVar= 1;
    [s_data_pos_filt, s_data_neg_filt, s_nf_pos_filt, s_nf_neg_filt, fs_data]= tdt_pink_helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);
    %     [nf_pos_data, nf_neg_data]= get_noisefloor_per_chin(data_dir, allfiles, tNF, fs_data);
    
    clean_s_data= nan(length(all_snrs)-1, 6);
    
    debugging= 0;
    if debugging
        nw_debug= 8;
        clf;
        ax(1)= subplot(121);
        hold on;
        plot_dpss_psd(s_nf_pos_filt, fs_data, 'nw', nw_debug);
        
        ax(2)= subplot(122);
        hold on;
        plot_dpss_psd(s_nf_neg_filt, fs_data, 'nw', nw_debug);
    end
    
    
    parfor snrVar= 2:length(all_snrs)
        curSNR= all_snrs(snrVar);
        
        [sn_data_pos_filt, sn_data_neg_filt, sn_nf_pos_filt, sn_nf_neg_filt, fs_data]= tdt_pink_helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);
        
        %         if ~debugging
        %             sn_nf_pos_filt= s_nf_pos_filt;
        %             sn_nf_neg_filt= s_nf_neg_filt;
        %         end
        
        figHan= 1;
        fName= strrep(sprintf('Q%d_nh_SNR%d_sn', chinID, curSNR), '-', 'm');
        ttlStr= sprintf('Q%d,NH,SNR %d', chinID, curSNR);
        
        if debugging
            subplot(121);
            plot_dpss_psd(sn_nf_pos_filt, fs_data, 'nw', nw_debug);
            subplot(122);
            plot_dpss_psd(sn_nf_neg_filt, fs_data, 'nw', nw_debug);
            
            if snrVar==length(all_snrs)
                legend(num2str(all_snrs'))
                linkaxes(ax);
            end
            
        else
            
            PSD_struct= create_panel_plot_s_vs_sn_tdt(figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, s_data_pos_filt, s_data_neg_filt, ...
                sn_nf_pos_filt, sn_nf_neg_filt, s_nf_pos_filt, s_nf_neg_filt, ...
                tStart, tEnd, fig_save_dir_subdir, fName, plotNF, saveFigs);
            
            raw_power_env(snrVar)   = PSD_struct.raw.ENV.SN;
            raw_power_tfs(snrVar)   = PSD_struct.raw.TFS.SN;
            
            frac_power_env(snrVar)  = PSD_struct.frac.ENV.SN;
            frac_power_tfs(snrVar)  = PSD_struct.frac.TFS.SN;
            
            raw_power_nf_env(snrVar)= PSD_struct.raw.ENV.NF_SN;
            raw_power_nf_tfs(snrVar)= PSD_struct.raw.TFS.NF_SN;
            
            raw_power_nf_env_CI(snrVar, :)= [PSD_struct.raw.ENV.NF_SN_CI_low PSD_struct.raw.ENV.NF_SN_CI_hi];
            raw_power_nf_tfs_CI(snrVar, :)= [PSD_struct.raw.TFS.NF_SN_CI_low PSD_struct.raw.TFS.NF_SN_CI_hi];
            
            clean_s_data(snrVar-1, :) = [PSD_struct.raw.ENV.S, PSD_struct.raw.TFS.S, PSD_struct.frac.ENV.S, PSD_struct.raw.TFS.S, PSD_struct.raw.ENV.NF_S, PSD_struct.raw.ENV.NF_S];
        end
        
    end
    clean_s_data= unique(clean_s_data, 'rows');
    
    raw_power_env(1)   = clean_s_data(1);
    raw_power_tfs(1)   = clean_s_data(2);
    frac_power_env(1)  = clean_s_data(3);
    frac_power_tfs(1)  = clean_s_data(4);
    
    raw_power_nf_env(1)= raw_power_nf_env(end-1);
    raw_power_nf_tfs(1)= raw_power_nf_tfs(end);
    
    
    raw_nf_power= -sqrt(raw_power_nf_env.*raw_power_nf_tfs);
    raw_power_nf_CI= -sqrt(raw_power_nf_env_CI.*raw_power_nf_tfs_CI);
    %%
    figure(11);
    hold on;
    mrkSize= 15;
    lw= 4;
    fSize= 24;
    [~, sort_inds]= sort(all_snrs);
    
    xlabel_str= cellfun(@(x) num2str(x), num2cell(all_snrs(sort_inds)), 'uniformoutput', false);
    xlabel_str= strrep(xlabel_str, '120', 'quiet');
    nf.min= nanmin(raw_power_nf_CI(:));
    nf.max= nanmax(raw_power_nf_CI(:));
    
    subplot(121);
    hold on;
    if plotNF
        l1(chinVar, 1)= plot(1:size(raw_power_env,1), raw_power_env(sort_inds), 'marker', mrkOrder(chinVar), ...
            'Color' , [.1 .1 1], 'markersize', mrkSize, 'linew', lw); % , 'DisplayName', sprintf('Q%d|ENV', chinID)
        l1(chinVar, 2)= plot(1:size(raw_power_tfs,1), raw_power_tfs(sort_inds), 'marker', mrkOrder(chinVar), ...
            'Color' , [1 .1 .1], 'markersize', mrkSize, 'linew', lw); % , 'DisplayName', sprintf('Q%d|TFS', chinID)
    else
        l1(chinVar, 1)= plot(1:size(raw_power_env,1), raw_power_env(sort_inds(1:end)), 'marker', mrkOrder(chinVar), ...
            'Color' , get_color('b'), 'markersize', mrkSize, 'linew', lw); % , 'DisplayName', sprintf('Q%d|ENV', chinID)
        l1(chinVar, 2)= plot(1:size(raw_power_tfs,1), raw_power_tfs(sort_inds(1:end)), 'marker', mrkOrder(chinVar), ...
            'Color' , get_color('r'), 'markersize', mrkSize, 'linew', lw); %, 'DisplayName', sprintf('Q%d|TFS', chinID)
    end
    
    subplot(122);
    hold on;
    if plotNF
        l1(chinVar, 1)= plot(1:size(raw_power_env,1), raw_power_env(sort_inds), 'marker', mrkOrder(chinVar), ...
            'Color' , [.1 .1 1], 'markersize', mrkSize, 'linew', lw); % , 'DisplayName', sprintf('Q%d|ENV', chinID)
        l1(chinVar, 2)= plot(1:size(raw_power_tfs,1), raw_power_tfs(sort_inds), 'marker', mrkOrder(chinVar), ...
            'Color' , [1 .1 .1], 'markersize', mrkSize, 'linew', lw); % , 'DisplayName', sprintf('Q%d|TFS', chinID)
    else
        l1(chinVar, 1)= plot(1:size(raw_power_env,1)-1, raw_power_env(sort_inds(1:end-1))-raw_power_env(sort_inds(end)), 'marker', mrkOrder(chinVar), ...
            'Color' , get_color('b'), 'markersize', mrkSize, 'linew', lw); % , 'DisplayName', sprintf('Q%d|ENV', chinID)
        l1(chinVar, 2)= plot(1:size(raw_power_tfs,1)-1, raw_power_tfs(sort_inds(1:end-1))-raw_power_tfs(sort_inds(end)), 'marker', mrkOrder(chinVar), ...
            'Color' , get_color('r'), 'markersize', mrkSize, 'linew', lw); %, 'DisplayName', sprintf('Q%d|TFS', chinID)
    end



    if plotNF
        warning('Not for NF - since considering drop');
        l1(chinVar, 3)= plot(1:size(raw_nf_power,1), raw_nf_power(sort_inds), 'marker', mrkOrder(chinVar), ...
            'Color' , [.1 1 .1], 'markersize', mrkSize, 'linew', lw); % , 'DisplayName', sprintf('Q%d|NF', chinID)
        
        % fl= fill([1:size(raw_power_env,1) fliplr(1:size(raw_power_env,1))], [raw_power_nf_CI(:,1)' fliplr(raw_power_nf_CI(:,2)')], 'g');
        fl= fill([1 size(raw_nf_power,1) size(raw_nf_power,1) 1 1], [nf.min nf.min nf.max nf.max nf.min], 'g');
        fl.FaceAlpha= .25;
        fl.EdgeAlpha= 0;
    end
    
    %     legend(l1, 'ENV', 'NF', 'Location', 'northwest');
    
    %     legend(l1, 'ENV', 'TFS', 'Location', 'northwest');
    grid on;
%     title(['NH- Q' num2str(chinID)]);
    
    
    
    if saveData
        save([data_save_dir 'SFRpink500_maskedPowers_Q' num2str(chinID) '.mat'], 'all_snrs', 'raw_power_nf_CI', ...
            'raw_power_env', 'frac_power_env', 'raw_power_tfs', 'frac_power_tfs', ...
            'raw_power_nf_env', 'raw_power_nf_tfs', 'raw_power_nf_env_CI', 'raw_power_nf_tfs_CI');
    end
end

tick_len= [.02 .02];

subplot(121);
set(gca,'fontsize', fSize, 'xtick', 1:size(raw_power_env,1), 'LineWidth', 1.5, 'TickLength', tick_len);
xticklabels(gca, xlabel_str(1:end));
ylabel('Total Power (dB)');
xlabel('HP masker SNR (dB)');
xlim([1 5])
grid off

subplot(122);
set(gca,'fontsize', fSize, 'xtick', 1:size(raw_power_env,1)-1, 'xticklabel', xlabel_str(1:end-1), 'LineWidth', 1.5, 'TickLength', tick_len);
ylabel('\DeltaPower re. quiet (dB)');
% xlabel('High-pass pink masker SNR (dB)');
grid off

lgHan(1)= plot(nan, nan, '-', 'color', get_color('r'), 'markersize', mrkSize, 'linew', lw);
lgHan(2)= plot(nan, nan, '-', 'color', get_color('b'), 'markersize', mrkSize, 'linew', lw);
lg = legend(lgHan, 'TFS', 'ENV');
lg.Location = 'southeast';
lg.Box= 'off';

add_subplot_letter(1, 2, fSize);
fName_summary= sprintf('NH_pink_summary');
set(gcf, 'units', 'centimeters', 'position', [40 3 30 16]);
% ylim([-55 -40])
if saveFinalFig
    saveas(gcf, [fig_save_dir_subdir fName_summary], 'png');
end
if saveLatex
    saveas(gcf, [LatexDir fName_summary], 'epsc');
end


