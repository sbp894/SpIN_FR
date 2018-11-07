%%
% Converts mfiles from NEL to MAT files


%%
function resample_FFR()

FsNew=16e3;

CodesDir=pwd;
addpath(CodesDir);


% NELDataRepository='R:\Users\Satya\SP\NELData\';
NELDataRepository='D:\Study Stuff\Matlab\NelData\MATData\';
MATDataRepository='D:\Study Stuff\Matlab\NelData\FFRMATDataResampled\';

DataDir=uigetdir(NELDataRepository);
OutDir=[MATDataRepository DataDir(length(fileparts(DataDir))+2:end) filesep];

mkdir(OutDir);

cd(DataDir);
allfiles=dir();

for file_var=1:length(allfiles)
    mfilename=allfiles(file_var).name;
    
    if strcmp(mfilename(1),'.') % Don't copy system dirs
        %root dirs
    elseif  ~isempty(strfind(mfilename, 'FFR')) % Resample
        load(mfilename);
        FsOld=data.Stimuli.RPsamprate_Hz;
        ADDataFields=fieldnames(data.AD_Data);
        
        for fieldVar=1:length(ADDataFields)
            if ~isempty(strfind(ADDataFields{fieldVar}, 'AD'))
                oldADdata=data.AD_Data.(ADDataFields{fieldVar});
                for repVar=1:length(oldADdata)
                    data.AD_Data.(ADDataFields{fieldVar}){repVar}=resample(oldADdata{repVar}, FsNew, round(FsOld));
                    data.Stimuli.FsNew=FsNew;
                end
            end
        end
        matfilename=[OutDir mfilename];
        save(matfilename,'data');
        
    elseif allfiles(file_var).isdir  % Copy directories
        copyfile(mfilename,[OutDir mfilename filesep]);
    else % Copy other files
        copyfile(mfilename,OutDir);
    end
end

cd(CodesDir);