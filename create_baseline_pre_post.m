clear;
clc;
close all;


saveFigs= 0;
allChins= [358 360 366 367 369 370];
% allChins= [358 360 366 367 370];

DataDir= '/media/parida/DATAPART1/Matlab/ABR/Output/';
figOutDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/';

LatexDir= '/home/parida/Dropbox/Seminars/SHRP_Feb19/figures/';
LatexDir_prelims= '/home/parida/Dropbox/Academics/Prelims/slides/figures/';

types= {'NH', 'HI'};
freqs2use= [.5 1 2 4 8]*1e3;

thresh_data.nh.z= {};
thresh_data.hi.z= {};
thresh_data.nh.amp= {};
thresh_data.hi.amp= {};

for chinVar= 1:length(allChins)
    cur_ChinID= allChins(chinVar);
    for typeVar= 1:length(types)
        cur_type= types{typeVar};
        curDir= dir(sprintf('%sQ%d_%s*', DataDir, cur_ChinID, cur_type));
        if ~isempty(curDir)
            fprintf('Using %s\n', curDir.name);
            curData_fName= dir(sprintf('%s%s%s*.mat', DataDir, curDir.name, filesep));
            curData= load(sprintf('%s%s%s%s', DataDir, curDir.name, filesep, curData_fName.name));
            curData= curData.abrs;
            [~, inds]= ismember(freqs2use, curData.thresholds(:,1));
            if strcmp(cur_type, 'NH')
                thresh_data.nh.z= [thresh_data.nh.z, curData.thresholds(inds,2)];
                %                 thresh_data.nh.amp= [thresh_data.nh.amp, curData.thresholds(:,3)];
            elseif strcmp(cur_type, 'HI')
                thresh_data.hi.z= [thresh_data.hi.z, curData.thresholds(inds,2)];
                %                 thresh_data.hi.amp= [thresh_data.hi.amp, curData.thresholds(:,3)];
            end
            
        else
            error('huh?');
        end
        
    end
end

%% plot
mrkSize= 20;
lw= 3;
fSize= 18;

figure(1);
clf;

subplot(211);
nh_data= cell2mat(thresh_data.nh.z)';
hi_data= cell2mat(thresh_data.hi.z)';
outlier_ind= 5;
reg_ind= setxor(1:length(allChins), outlier_ind);

plot(freqs2use, nh_data(reg_ind, :), '-bd', freqs2use, nh_data(outlier_ind, :), '-cd','markersize', mrkSize, 'linew', lw);
hold on;
plot(freqs2use, hi_data(reg_ind, :), '-ro', freqs2use, hi_data(outlier_ind, :), '-mo','markersize', mrkSize, 'linew', lw);
grid on;
set(gca, 'xscale', 'log', 'xtick', freqs2use);
xlim([.4 10]*1e3);
xlabel('Frequency (Hz)');
ylabel('ABR threshold (dB)');
set(gca, 'fontsize', fSize);

lg= plot(nan, nan, '-bd', nan, nan, '-ro','markersize', mrkSize, 'linew', lw);
legend(lg, 'NH', 'HI');


subplot(223);
plot(freqs2use, hi_data(reg_ind, :)-nh_data(reg_ind, :), '-bd','markersize', mrkSize, 'linew', lw);
hold on;
plot(freqs2use, hi_data(outlier_ind, :)-nh_data(outlier_ind, :), '-co','markersize', mrkSize, 'linew', lw);
grid on;
set(gca, 'xscale', 'log', 'xtick', freqs2use);
xlim([.4 10]*1e3);
xlabel('Frequency (Hz)');
ylabel('ABR threshold shift (dB)');
set(gca, 'fontsize', fSize);
lg2= plot(nan, nan, '-co','markersize', mrkSize, 'linew', lw);
legend(lg2, 'Q369', 'location', 'southwest');


subplot(224);
plot(reg_ind, (hi_data(reg_ind, :)-nh_data(reg_ind, :)), 'bd', 'markersize', mrkSize, 'linew', lw);
hold on;
plot(outlier_ind, (hi_data(outlier_ind, :)-nh_data(outlier_ind, :)), 'co','markersize', mrkSize, 'linew', lw);
grid on;
set(gca, 'xtick', 1:length(allChins), 'xticklabel', cellfun(@(x) ['Q' num2str(x)], num2cell(allChins), 'uniformoutput', false));
xlim([0.5 length(allChins)+0.5]);
xlabel('Animal ID');
ylabel('ABR threshold shift (dB)');
set(gca, 'fontsize', fSize);

set(gcf, 'units', 'normalized', 'position', [.1 .1 .8 .8]);
fName= 'baselines_for_sfr';
% saveas(gcf, [figOutDir fName], 'png');

% saveas(gcf, [LatexDir fName], 'epsc');

%%
figure(2);
clf;

nh_data= cell2mat(thresh_data.nh.z)';
hi_data= cell2mat(thresh_data.hi.z)';
outlier_ind= 5;
reg_ind= setxor(1:length(allChins), outlier_ind);

plot(freqs2use, nh_data, '-bd','markersize', mrkSize, 'linew', lw);
hold on;
plot(freqs2use, hi_data, '-ro','markersize', mrkSize, 'linew', lw);
grid on;
set(gca, 'xscale', 'log', 'xtick', freqs2use);
xlim([.4 10]*1e3);
xlabel('Frequency (Hz)');
ylabel('ABR threshold (dB)');
set(gca, 'fontsize', fSize);

lg= plot(nan, nan, '-bd', nan, nan, '-ro','markersize', mrkSize, 'linew', lw);
legend(lg, 'NH', 'HI', 'Location', 'southeast');

ylim([-2 50]);
set(gcf, 'Units', 'inches', 'Position', [1 1 6 4]);
fName2= 'abr_thresh_ffr';
if saveFigs
    saveas(gcf, [LatexDir_prelims fName2], 'epsc');
end