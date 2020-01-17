clear;
clc;

CodesDirs= {'/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/helper', ...
    '/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/continuous'};
addpath(CodesDirs{:});

plotVar=1;
freq_output_spread=0;

outFigDir= '/home/parida/Dropbox/Articles/Loss_of_tonotopy_in_HI_FFR/figures/';
saveFig= 1;
fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/');
if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end

[stim, fsOld]=audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2017_09_09-Q321_AN_NH/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_N.wav');
fsStim= 10e3;
stim= gen_resample(stim, fsOld, fsStim);

tStim= (1:length(stim))/fsStim;

nh.Dir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2019_03_19-Q371_SFR_artifact_NH/';
nh.allfiles= dir([nh.Dir '*.mat']);

nh.art_file= nh.allfiles(contains({nh.allfiles.name}, '_artifact')).name;
nh.reg_file= nh.allfiles(contains({nh.allfiles.name}, 'snr_120') & ~contains({nh.allfiles.name}, '_artifact')).name;

temp= load([nh.Dir nh.art_file]);
fs_old= temp.data.Stimuli.RPsamprate_Hz;
fs_data= 5e3;

temp.data.AD_Data.AD_Avg_PO_V{1}= gen_resample(temp.data.AD_Data.AD_Avg_PO_V{1}, fs_old, fs_data);
temp.data.AD_Data.AD_Avg_NP_V{1}= gen_resample(temp.data.AD_Data.AD_Avg_NP_V{1}, fs_old, fs_data);

fMax= fs_data/2;

clf;
nh.art_data.pos= remove_artifact_ffr(temp.data.AD_Data.AD_Avg_PO_V{1}, fs_data, plotVar, fMax, freq_output_spread);
clf;
nh.art_data.neg= remove_artifact_ffr(temp.data.AD_Data.AD_Avg_NP_V{1}, fs_data, plotVar, fMax, freq_output_spread);
clf;

temp= load([nh.Dir nh.reg_file]);
fs_old= temp.data.Stimuli.RPsamprate_Hz;
temp.data.AD_Data.AD_Avg_PO_V{1}= gen_resample(temp.data.AD_Data.AD_Avg_PO_V{1}, fs_old, fs_data);
temp.data.AD_Data.AD_Avg_NP_V{1}= gen_resample(temp.data.AD_Data.AD_Avg_NP_V{1}, fs_old, fs_data);


nh.reg_data.pos= remove_artifact_ffr(temp.data.AD_Data.AD_Avg_PO_V{1}, fs_data, plotVar, fMax, freq_output_spread);
clf;
nh.reg_data.neg= remove_artifact_ffr(temp.data.AD_Data.AD_Avg_NP_V{1}, fs_data, plotVar, fMax, freq_output_spread);
clf;

t_data= (1:length(nh.reg_data.pos))/fs_data;

HalfPowerFrequency1= 70;
HalfPowerFrequency2= 2e3;
N_bp_half= 4;

bpFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);

%
nh.reg_data.pos= filtfilt(bpFilt, nh.reg_data.pos);
nh.reg_data.neg= filtfilt(bpFilt, nh.reg_data.neg);
nh.reg_data.env= (nh.reg_data.pos+nh.reg_data.neg)/2;
nh.reg_data.tfs= (nh.reg_data.pos-nh.reg_data.neg)/2;

%
nh.art_data.pos= filtfilt(bpFilt, nh.art_data.pos);
nh.art_data.neg= filtfilt(bpFilt, nh.art_data.neg);
nh.art_data.env= (nh.art_data.pos+nh.art_data.neg)/2;
nh.art_data.tfs= (nh.art_data.pos-nh.art_data.neg)/2;


%%
figure(1);
clf;


fSize= 24;
lg_fSize= 20;
gain = 20e3/10;
tickVals= [0 .5 1 1.5];

stim= .5*stim/rms(stim)*rms(nh.reg_data.pos/gain*1e6);

Amax=25;
subplot(121);
set_colblind_order();


hold on;
% grid on;
tHan(1)= plot(t_data, nh.reg_data.pos/gain*1e6);
tHan(2)= plot(t_data, -Amax+nh.art_data.pos/gain*1e6);
% set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1)
plot(tStim, -2*Amax+stim, 'k');
% ylabel({'{\bfFFR +VE}'; 'Amplitude (arb)'});
ylabel('Amplitude (\muV)');
% title('Amplitude');
% xlabel('time (s)');
tick_len= [.025 .025];
set(gca, 'FontSize', fSize, 'LineWidth', 1.5, 'TickLength', tick_len, 'XTick', tickVals)
xlabel('Time (sec)');


% subplot(223);
% hold on;
% grid on;
% plot(t_data, nh.reg_data.neg);
% plot(t_data, -Amax+nh.art_data.neg');
% % set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')+1)
% plot(tStim, -2*Amax-.1*stim, 'k');
% ylabel({'{\bfFFR -VE}'; 'Amplitude (arb)'});
% set(gca, 'FontSize', fSize, 'YTick', [])
% % set(gca, 'FontSize', fSize)
l0= plot(nan, nan, 'k', 'linew', 2);

% subplot(222);
% hold on;
% grid on;
% plot(t_data, nh.reg_data.env);
% plot(t_data, nh.reg_data.tfs);
%
% plot(t_data, -Amax+nh.art_data.env);
% plot(t_data, -Amax+nh.art_data.tfs);


yRange= 42;
nw=5;
nfft= round(fs_data);
xtick_vals= [20 100 500 2e3];
xtick_labs= cellfun(@(x) num2str(x), num2cell(xtick_vals), 'UniformOutput', false);

% title('A', 'horizontalAlignment', 'left');

stim= stim/400;

subplot(122);
set_colblind_order();

hold on;
[~,~,l3]= plot_dpss_psd(stim, fsStim, 'nfft', nfft, 'nw', nw, 'plotconf', false);
set(gca, 'ColorOrderIndex', 1);
[~,~,l1]= plot_dpss_psd(nh.reg_data.pos, fs_data, 'nfft', nfft, 'nw', nw, 'plotconf', true);
[~,~,l2]= plot_dpss_psd(nh.art_data.pos, fs_data, 'nfft', nfft, 'nw', nw, 'plotconf', true);
set(l3, 'color', 'k', 'LineStyle', '-', 'LineWidth', 2);
axis tight;
yl= ylim;
ylimHard= [max(yl)-yRange+3 max(yl)+3];
ylim(ylimHard);
xlim([20 2e3]);
% legend([l1(1) l2(1)], 'With eartip', 'Without eartip')
set(gca, 'FontSize', fSize, 'XTick', xtick_vals, 'XTickLabel', xtick_labs, 'LineWidth', 1.5, 'TickLength', tick_len)
ylabel('PSD (dB/Hz)');
title('')
xlabel('')
grid off;
% xlabel('Frequency (Hz)')
% ttlHan= title('B');
% set(ttlHan, 'horizontalAlignment', 'left');

xVert= [90, 500];

matGreen= get_color('g');
% plot(repmat(xVert, 2, 1), repmat(ylimHard', 1, 2), 'Color', matGreen, 'LineWidth', 3, 'LineStyle', '-');
plot(xVert, [min(ylimHard) min(ylimHard)]+2, 'Color', matGreen, 'LineWidth', 10, 'LineStyle', '-');

lg=legend([l1(1) l2(1) l0], '{\bf+}eartip', '{\bf-}eartip', 'stim', 'Location', 'northeast');
lg.FontSize=lg_fSize;
lg.Box= 'off';

%
% subplot(224);
% hold on;
% plot_dpss_psd(nh.reg_data.env, fs_data, 'nfft', nfft, 'nw', nw);
% plot_dpss_psd(nh.reg_data.tfs, fs_data, 'nfft', nfft, 'nw', nw);
% [~,~,l3]= plot_dpss_psd(nh.art_data.env, fs_data, 'nfft', nfft, 'nw', nw);
% set(l3, 'color', 'g', 'LineStyle', '--');
% [~,~,l4]= plot_dpss_psd(nh.art_data.tfs, fs_data, 'nfft', nfft, 'nw', nw);
% legend('data-env', 'data-tfs', 'art-env', 'art-tfs')
% set(l4, 'color', 'c', 'LineStyle', '--');
% axis tight;
% yl= ylim;
% ylim([max(yl)-yRange max(yl)]);

% subplot(224);
% hold on;
% [~,~,~]=plot_dpss_psd(nh.reg_data.neg, fs_data, 'nfft', nfft, 'nw', nw, 'plotconf', true);
% [~,~,~]= plot_dpss_psd(nh.art_data.neg, fs_data, 'nfft', nfft, 'nw', nw, 'plotconf', true);
% % set(l2(2), 'color', 'g', 'LineStyle', '-', 'LineWidth', 1.2);
% axis tight;
% yl= ylim;
% ylim([max(yl)-yRange max(yl)]);
% xlim([20 2e3]);
% set(gca, 'FontSize', fSize, 'XTick', xtick_vals, 'XTickLabel', xtick_labs)
% title('');
% ylabel('PSD (dB/Hz)');
xlabel('Frequency (Hz)')

set(gcf, 'Units', 'inches', 'Position', [1 1 11 4.5]);
% saveas(gcf, [fig_save_dir 'compare_no_eartip'], 'png')

txtHan= add_subplot_letter(1, 2, 30);

if saveFig
    saveas(gcf, [outFigDir 'compare_no_eartip'], 'epsc')
end

rmpath(CodesDirs{:});