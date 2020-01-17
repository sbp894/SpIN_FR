clear;
clc;

chinID= 379;

if chinID==373
    picNums_to_combine= [2 7];
elseif chinID==374
    picNums_to_combine= [4 7];
elseif chinID==379
    picNums_to_combine= [2 7];
end



RootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
data_dir= dir([RootDataDir '*Q' num2str(chinID) '*pink500*']);
data_dir= [RootDataDir data_dir.name filesep];



all_p_files= dir([data_dir 'p*SFR*.mat']);

picNums= cell2mat(cellfun(@(x) sscanf(x, 'p%04d*'), {all_p_files.name}, 'UniformOutput', false));

tempData= cell(length(picNums_to_combine), 1);
nPairs= nan(length(picNums_to_combine), 1);
nPairs_actual= nan(length(picNums_to_combine), 1);

for fileVar= 1:length(picNums_to_combine)
    curPicNum= picNums_to_combine(fileVar);
    curFile= all_p_files(picNums==curPicNum).name;
    xx= load([data_dir curFile]);
    xx= xx.data;
    
    fprintf('Using %s \n', curFile);
    
    tempData{fileVar}=  xx.AD_Data.AD_All_V;
    nPairs(fileVar)= xx.Stimuli.RunLevels_params.nPairs;
    nPairs_actual(fileVar) = xx.Stimuli.RunLevels_params.nPairs_actual;
    
    if fileVar==1
        data= xx;
    else 
        data.AD_Data.AD_All_V= [data.AD_Data.AD_All_V, xx.AD_Data.AD_All_V];
        data.Stimuli.RunLevels_params.nPairs= data.Stimuli.RunLevels_params.nPairs+xx.Stimuli.RunLevels_params.nPairs;
        data.Stimuli.RunLevels_params.nPairs_actual= data.Stimuli.RunLevels_params.nPairs_actual+xx.Stimuli.RunLevels_params.nPairs_actual;
    end
end

data.AD_Data.AD_Avg_PO_V{1}= nanmean(cell2mat(data.AD_Data.AD_All_V(1:2:end)'), 1);
data.AD_Data.AD_Avg_NP_V{1}= nanmean(cell2mat(data.AD_Data.AD_All_V(2:2:end)'), 1);

max_PicNum= dir([data_dir 'p*']);
max_PicNum= max(cell2mat(cellfun(@(x) sscanf(x, 'p%04d*'), {max_PicNum.name}, 'UniformOutput', false)));
p_fName2Save= sprintf('%sp%04d_FFR_SNRenvSSN_Stim_S_P_atn10.mat', data_dir, max_PicNum+1);

save(p_fName2Save, 'data');

data.AD_Data= rmfield(data.AD_Data, 'AD_All_V');
a_fName2Save= sprintf('%sa%04d_FFR_SNRenvSSN_Stim_S_P_atn10.mat', data_dir, max_PicNum+1);
save(a_fName2Save, 'data');


warning('Check if filename is correct');