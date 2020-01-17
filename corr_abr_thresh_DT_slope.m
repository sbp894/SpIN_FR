clear;
clc;

ABR_DataDir= '/media/parida/DATAPART1/Matlab/ABR/Output/';

saveFig= 0;

fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/moving_segment_analysis_vowel/');
if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end

allChinData = load('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Summary_slopes/all_chins_data.mat');
allChinData = allChinData.allChinData;

freqs2use= [.5 1 2 4 8 0]*1e3;
midFreqs= [1 2 4]*1e3;

for iterVar= 1:length(allChinData)
    chinID= allChinData(iterVar).chinID;
    chinType= allChinData(iterVar).group;
    
    lf_power_audio= dbspl(sqrt(allChinData(iterVar).lf_power_audio));
    hf_power_audio= dbspl(sqrt(allChinData(iterVar).hf_power_audio));
    
    tfs_power_ffr= db(sqrt(allChinData(iterVar).tfs_power_ffr));
    env_power_ffr= db(sqrt(allChinData(iterVar).env_power_ffr));
    
    
    mdl_raw= fitlm(lf_power_audio, tfs_power_ffr);
    allChinData(iterVar).DTslope_raw = mdl_raw.Coefficients.Estimate(2);
    
    mdl_norm= fitlm(hf_power_audio - lf_power_audio, env_power_ffr - tfs_power_ffr);
    allChinData(iterVar).DTslope_norm = mdl_norm.Coefficients.Estimate(2);
    
    
    if any(strcmpi(chinType, {'PTS', 'HI'}))
        chinType= 'HI';
    end
    
    curDir= dir(sprintf('%sQ%d_%s*', ABR_DataDir, chinID, chinType));
    
    if numel(curDir)>1
        warning('Using the directory %s\n', curDir(1).name);
        curDir= curDir(1);
    end
    if ~isempty(curDir)
        
        fprintf('Using %s\n', curDir.name);
        
        curData_fName= dir(sprintf('%s%s%s*.mat', ABR_DataDir, curDir.name, filesep));
        curData= load(sprintf('%s%s%s%s', ABR_DataDir, curDir.name, filesep, curData_fName.name));
        curData= curData.abrs;
        [~, inds]= ismember(freqs2use, curData.thresholds(:,1));
        
        allChinData(iterVar).abr_thresh= nan(size(freqs2use));
        allChinData(iterVar).abr_thresh= curData.thresholds(find(inds),2); %#ok<FNDSB>
        
        allChinData(iterVar).midFreq_thresh= mean(allChinData(iterVar).abr_thresh(ismember(freqs2use, midFreqs)));
        
    else
        fprintf('%s \n ', '---');
        allChinData(iterVar).midFreq_thresh = nan;
    end
    
end

abr_thresh= [allChinData.midFreq_thresh];
dt_slope_raw= [allChinData.DTslope_raw]
dt_slope_norm= [allChinData.DTslope_norm]

%%
figure(2);
clf;
special_cond= [allChinData.chinID] == 369 & strcmp({allChinData.group}, 'PTS');

subplot(211)
hold on;
plot(abr_thresh, dt_slope_raw, 'd', 'linew', 2, 'markersize', 14)
plot(abr_thresh(special_cond), dt_slope_raw(special_cond), 'cd', 'linew', 3, 'markersize', 14)
ylabel('TFS vs LF');
title('Slope-Raw')
grid off;

subplot(212)
hold on;
plot(abr_thresh, dt_slope_norm, 'd', 'linew', 2, 'markersize', 14)
plot(abr_thresh(special_cond), dt_slope_norm(special_cond), 'cd', 'linew', 3, 'markersize', 14)

xlabel('mean ABR thresh (1,2,4k)');
ylabel('TFS/E vs LF/HF');
title('Slope-Normalized ');
grid off;

set(findall(gcf,'-property','FontSize'),'FontSize',16)

set(gcf, 'Units', 'inches', 'Position', [1 1 8 6]);

if saveFig
    saveas(gcf, [fig_save_dir 'thresh_vs_slope'], 'png');
end