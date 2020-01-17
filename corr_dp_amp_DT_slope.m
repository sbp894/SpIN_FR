clear;
clc;

ABR_DataDir= '/media/parida/DATAPART1/Matlab/ABR/Output/';
DPoae_rootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/Baselines/';
LatexDir= '/home/parida/Dropbox/Articles/Loss_of_tonotopy_in_HI_FFR/figures/';

CodesDir= '/media/parida/DATAPART1/Matlab/Screening';
addpath(CodesDir);

fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/moving_segment_analysis_vowel/');
if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end

allChinData = load('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Summary_slopes/all_chins_data.mat');
allChinData = allChinData.allChinData;

freqs2use= [.5 1 2 4 8 0]*1e3;
midFreqs= [1 2 4 8]*1e3;
fSize= 20;
saveFig= 1;

for iterVar= 1:length(allChinData)
    chinID= allChinData(iterVar).chinID;
    chinType= allChinData(iterVar).group;
    
    lf_power_audio= dbspl(sqrt(allChinData(iterVar).lf_power_audio));
    hf_power_audio= dbspl(sqrt(allChinData(iterVar).hf_power_audio));
    
    tfs_power_ffr= db(sqrt(allChinData(iterVar).tfs_power_ffr));
    env_power_ffr= db(sqrt(allChinData(iterVar).env_power_ffr));
    
    notNANinds= ~(isnan(lf_power_audio) | isnan(hf_power_audio) | isnan(tfs_power_ffr) | isnan(env_power_ffr));
    
    if any(notNANinds)
        
        mdl_raw= fitlm(lf_power_audio(notNANinds), tfs_power_ffr(notNANinds));
        allChinData(iterVar).DTslope_raw = mdl_raw.Coefficients.Estimate(2);
        
        mdl_norm= fitlm(hf_power_audio(notNANinds)- lf_power_audio(notNANinds), env_power_ffr(notNANinds) - tfs_power_ffr(notNANinds));
        allChinData(iterVar).DTslope_norm = mdl_norm.Coefficients.Estimate(2);
        
        if any(strcmpi(chinType, {'PTS', 'HI'}))
            chinType= 'HI';
        end
        
        %% Do ABR analyis
        curDir= dir(sprintf('%sQ%d_%s*', ABR_DataDir, chinID, chinType));
        if ~isempty(curDir)
            
            if numel(curDir)>1
                warning('Check that the right directory is selected');
                curDir= curDir(1);
            end
            
            fprintf('Using %s\n', curDir.name);
            
            curData_fName= dir(sprintf('%s%s%s*.mat', ABR_DataDir, curDir.name, filesep));
            curData= load(sprintf('%s%s%s%s', ABR_DataDir, curDir.name, filesep, curData_fName.name));
            curData= curData.abrs;
            [~, inds]= ismember(freqs2use, curData.thresholds(:,1));
            
            allChinData(iterVar).abr_thresh= nan(size(freqs2use));
            allChinData(iterVar).abr_thresh= curData.thresholds(find(inds),2); %#ok<FNDSB>
            
            allChinData(iterVar).midFreq_thresh= nanmean(allChinData(iterVar).abr_thresh(ismember(freqs2use, midFreqs)));
            
        else
            fprintf('%s \n ', '---');
            allChinData(iterVar).midFreq_thresh = nan;
        end
        
        %% Do DPOAE analysis

        allChinDirs= dir([DPoae_rootDataDir '*' num2str(chinID) '*']);
        if~isempty(allChinDirs)
            if strcmp(chinType, 'NH')
                dirNum= find(contains({allChinDirs.name}, {'nh', 'pre'},'IgnoreCase',true));
            elseif strcmp(chinType, 'PTS') | strcmp(chinType, 'HI')
                dirNum= find(contains({allChinDirs.name}, {'pts', 'post', 'follow', 'hi'}));
            else
                dirNum= nan;
            end
            if ~isempty(dirNum) && any(~isnan(dirNum))
                if numel(dirNum)>1
                    warning('Check that the right directory is selected');
                    dirNum= dirNum(1);
                end
                DataDir= allChinDirs(dirNum).name;
                dpFile= dir([DPoae_rootDataDir DataDir filesep '*dpoae*']);
                dpFile= [DPoae_rootDataDir DataDir filesep dpFile(1).name];
                calibFile= get_lower_calibFile(dpFile);
                
                run(calibFile);
                calibData=ans;
                calibData=calibData.CalibData;
                
                out_DPOAE_data= my_dpoae_analysis(dpFile);
                dpData=[out_DPOAE_data.dp_amp];
                dp_freqs= [out_DPOAE_data.freq2];
                
                freq_range2consider= [min(midFreqs) max(midFreqs)];
                inds2consider= dp_freqs>freq_range2consider(1) & dp_freqs<freq_range2consider(2);
                allChinData(iterVar).midFreq_DPamp= nanmean(dpData(inds2consider));
            else
                allChinData(iterVar).midFreq_DPamp= nan;
            end
        else
            allChinData(iterVar).midFreq_DPamp= nan;
        end
    else
        allChinData(iterVar).midFreq_thresh= nan;
        allChinData(iterVar).midFreq_DPamp= nan;
        allChinData(iterVar).DTslope_raw= nan;
        allChinData(iterVar).DTslope_norm= nan;
    end
    
end

abr_thresh_db= [allChinData.midFreq_thresh];
dp_amp_db= [allChinData.midFreq_DPamp];
dt_slope_raw= [allChinData.DTslope_raw];
dt_slope_norm= [allChinData.DTslope_norm];

%%
plotYes= 1;
verbose= 1;
lw=3;
lw2= 4;

figure(2);
clf;
special_cond= [allChinData.chinID] == 369 & strcmp({allChinData.group}, 'PTS');

% subplot(221)
% hold on;
% plot(abr_thresh_db, dt_slope_raw, 'd', 'linew', lw, 'markersize', 14)
% plot(abr_thresh_db(special_cond), dt_slope_raw(special_cond), 'cd', 'linew', lw2, 'markersize', 14);
% plot_fitlm(abr_thresh_db, dt_slope_raw, plotYes, verbose)
% ylabel('TFS vs LF');
% title('ABR')
% grid off;
%
% subplot(222)
% hold on;
% plot(dp_amp_db, dt_slope_raw, 'd', 'linew', lw, 'markersize', 14)
% plot(dp_amp_db(special_cond), dt_slope_raw(special_cond), 'cd', 'linew', lw2, 'markersize', 14)
% plot_fitlm(dp_amp_db, dt_slope_raw, plotYes, verbose)
% title('DPOAE')
% grid off;

tick_len= [.025 .025];
subplot(121)
hold on;
plot(abr_thresh_db, dt_slope_norm, 'kd', 'linew', lw, 'markersize', 14)
plot(abr_thresh_db(special_cond), dt_slope_norm(special_cond), 'd', 'color', get_color('lr'), 'linew', lw2, 'markersize', 14)
plot_fitlm(abr_thresh_db, dt_slope_norm, plotYes, verbose)
% xlabel(sprintf('mean_{%.0f-%.0fkHz} ABR thresh (dB)', min(midFreqs)/1e3, max(midFreqs)/1e3));
xlabel('ABR Threshold (dB SPL)');
ylabel('Normalized DT slope');
set(gca, 'LineWidth', 1.5, 'TickLength', tick_len);
grid off;

subplot(122)
hold on;
plot(dp_amp_db, dt_slope_norm, 'kd', 'linew', lw, 'markersize', 14)
plot(dp_amp_db(special_cond), dt_slope_norm(special_cond), 'd', 'color', get_color('lr'), 'linew', lw2, 'markersize', 14)
plot_fitlm(dp_amp_db, dt_slope_norm, plotYes, verbose)
% xlabel(sprintf('mean_{%.0f-%.0fkHz} DP Amp (dB)', min(midFreqs)/1e3, max(midFreqs)/1e3));
xlabel('DP Amplitude (dB SPL) ');
grid off;
set(gca, 'LineWidth', 1.5, 'TickLength', tick_len);

set(findall(gcf,'-property','FontSize'),'FontSize',fSize)

set(gcf, 'Units', 'inches', 'Position', [1 1 11 5]);
add_subplot_letter(1, 2, 30)

if saveFig
    saveas(gcf, [fig_save_dir 'baseline_vs_slope'], 'png');
    saveas(gcf, [LatexDir 'baseline_vs_slope'], 'epsc');
end

rmpath(CodesDir);