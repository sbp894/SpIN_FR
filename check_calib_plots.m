clear;
clc;
clf;

AllChinIDs= [369 369 370 358 360 366 367];

tStart = 0; tEnd = 1.3;

fig_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/Calib/';
if ~isdir(fig_save_dir)
    mkdir(fig_save_dir);
end

xtickVals= [.2 .5 1 2 4 7 10];
fSize=16;
lw=2;
saveFigs= 1;
ax= nan(length(AllChinIDs), 1);
legStr= cell(length(AllChinIDs), 1);

for chinVar= 1:length(AllChinIDs)
    chinID= AllChinIDs(chinVar);
    RootDataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
    
    allFiles= dir([RootDataDir '*' num2str(chinID) '*SFR*']);
    
    if isempty(allFiles)
        error('No dir. what to do?');
    elseif length(allFiles)>1
        warning('there are multiple dirs. choosing last one');
        warning('bad!!! ');
        if chinVar==1
            data_dir= [RootDataDir allFiles(1).name filesep];
        else
            data_dir= [RootDataDir allFiles(end).name filesep];
        end
    else
        data_dir= [RootDataDir allFiles.name filesep];
    end
    
    [~, legStr{chinVar}]= fileparts(data_dir(1:end-1));
    
    files_in_dir= dir([data_dir '*calib*.mat']);
    if length(files_in_dir)==1
        temp= load([data_dir  files_in_dir.name]);
    elseif length(files_in_dir)>1
        temp= load([data_dir  files_in_dir(end).name]);
        warning('using last calib file');
    end
    calibData= temp.data.CalibData;
    figure(1);
    hold on;
    %     clf;
    ax(chinVar)= semilogx(calibData(:,1), calibData(:,2), 'linew', lw);
end

legStr= cellfun(@(x) strrep(x, '_', '-'), legStr, 'uniformoutput', false);
legStr= cellfun(@(x) x((find(ismember(legStr{1}, 'Q'))):end), legStr, 'uniformoutput', false);
set(gca, 'xtick', xtickVals, 'fontsize', fSize, 'xscale', 'log');
axis tight;
grid on;
xlabel('frequency (kHz)');
ylabel('Calib data');
xlim([0.2 10]);
title('Calib plot');

lg= legend(ax, legStr);
set(gcf, 'units', 'normalized', 'position', [0 0 1 1]);

fName= 'calib_allChins';
saveas(gcf, [fig_save_dir fName], 'tiff');