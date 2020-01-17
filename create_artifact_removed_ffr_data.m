clear;
clc;

chinID= 377;
verbose= 1;

rng(0);

CodesDirs= {'/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/helper', ...
    '/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/continuous'};
addpath(CodesDirs{:});

RootData_load_dir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
RootData_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Artifact_Removed_FFR/';
freqMax= 3e3;

if ~isfolder(RootData_save_dir)
    mkdir(RootData_save_dir);
end


search_data_dir= [dir([RootData_load_dir '*Q' num2str(chinID) '*FFR*']); dir([RootData_load_dir '*Q' num2str(chinID) '*SFR*'])];

if length(search_data_dir)>1
    for fileVar=1:length(search_data_dir)
        fprintf('(%d) %s \n', fileVar, search_data_dir(fileVar).name);
    end
    dirNum= input('Which directory?');
else
    dirNum= 1;
end


load_data_dir= [RootData_load_dir search_data_dir(dirNum).name filesep];
save_data_dir= [RootData_save_dir search_data_dir(dirNum).name filesep];

if ~isfolder(save_data_dir)
    mkdir(save_data_dir);
    flagProceed= 1;
else
    userIn= questdlg('Skip or overwrite?', 'Directory already exists.', 'Skip', 'Overwrite', 'Skip');
    if strcmp(userIn, 'Skip')
        flagProceed= 0;
    elseif strcmp(userIn, 'Overwrite')
        flagProceed= 1;
    end
end

if flagProceed
    fprintf('Choosing -------------------\n loading dir:%s \n saving dir: %s \n', load_data_dir, save_data_dir);
    fprintf('--------------------------------------------------------\n--------------------------------------------------------\n');
    
    all_afiles= [dir([load_data_dir 'p*FFR*.mat*']); dir([load_data_dir 'p*SFR*.mat*'])]; %[dir([load_data_dir 'a*FFR*.mat*']); dir([load_data_dir 'a*SFR*.mat*'])];
    
    plotVar= 0;
    
    for fileVar= 1:length(all_afiles)
        fprintf('[%d/%d] .............. Working on %s \n', fileVar, length(all_afiles), all_afiles(fileVar).name);
        
        cur_p_fName= [load_data_dir all_afiles(fileVar).name];
        new_p_fName= [save_data_dir all_afiles(fileVar).name];
        new_a_fName= [save_data_dir 'a' all_afiles(fileVar).name(2:end)];
        
        if ~exist(new_p_fName, 'file')
            data= load(cur_p_fName);
            data=data.data;
            
            if isfield(data.AD_Data, 'AD_All_V')
                temp_pos= cell2mat(data.AD_Data.AD_All_V(1:2:end)');
                temp_neg= cell2mat(data.AD_Data.AD_All_V(2:2:end)');
                rand_inds_pos= randsample(size(temp_pos,1), round(size(temp_pos,1)/2));
                rand_inds_neg= setxor(1:size(temp_pos,1), rand_inds_pos);
                
                data.AD_Data.AD_NF_PO_V{1}= nanmean(temp_pos(rand_inds_pos, :), 1) - nanmean(temp_pos(rand_inds_neg, :), 1);
                data.AD_Data.AD_NF_NP_V{1}= nanmean(temp_neg(rand_inds_pos, :), 1) - nanmean(temp_neg(rand_inds_neg, :), 1);
                
                data.AD_Data.AD_NF_PO_V{1} = remove_artifact_ffr(data.AD_Data.AD_NF_PO_V{1}, data.Stimuli.RPsamprate_Hz, plotVar, freqMax);
                data.AD_Data.AD_NF_NP_V{1} = remove_artifact_ffr(data.AD_Data.AD_NF_NP_V{1}, data.Stimuli.RPsamprate_Hz, plotVar, freqMax);
            end
            
            data.AD_Data.AD_Avg_PO_V{1}= remove_artifact_ffr(data.AD_Data.AD_Avg_PO_V{1}, data.Stimuli.RPsamprate_Hz, plotVar, freqMax);
            data.AD_Data.AD_Avg_NP_V{1}= remove_artifact_ffr(data.AD_Data.AD_Avg_NP_V{1}, data.Stimuli.RPsamprate_Hz, plotVar, freqMax);
            
            save(new_p_fName, 'data');
            
            data.AD_Data= rmfield(data.AD_Data, 'AD_All_V');
            save(new_a_fName, 'data');
        end
    end
    fprintf('------------------Done--------------------\n');
end
rmpath(CodesDirs{:});