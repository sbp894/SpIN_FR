function PSD_struct= create_panel_plot_s_vs_sn...
    (figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, s_data_pos_filt, s_data_neg_filt, tStart, tEnd, fig_save_dir, fName, ttlStr, saveFigs)

sn_data_env= (sn_data_pos_filt+sn_data_neg_filt)/2;
sn_data_tfs= (sn_data_pos_filt-sn_data_neg_filt)/2;

s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;

if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end

t_sig= (1:length(sig))/fs_sig;
t_data= (1:length(sn_data_pos_filt))/fs_data;

inds2use_stim= t_sig>tStart & t_sig<tEnd;
inds2use_data= t_data>tStart & t_data<tEnd;
inds2use_nfloor= t_data>length(sig)/fs_sig;

yRange_stim= 50;
yRange_ffr= 20;
lw2= 2;
nw=5;
nSProws=1;
nSPcols=2;
fSize= 20;

figure(figHan);
clf;
co= get(gca, 'colororder');
lw3= 3;
ylHard= [-50-yRange_stim -50];
xlHard= [75 1.5e3];
stimColor= .4*[1 1 1];

%%
xtick_vals= [10 100 1e3 4e3];
xtick_labs= cellfun(@(x) num2str(x), num2cell(xtick_vals), 'UniformOutput', false);
subplot(nSProws, nSPcols, 1);
yyaxis left;
[Pxx, ~, ax]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
ylim([max(Pxx)+5-yRange_stim max(Pxx)+5]);
ax_left(1)= gca;
ylabel('Stim-PSD (dB/Hz)');
set(ax, 'color', stimColor);
set(gca,'ycolor', stimColor);

yyaxis right;
hold on;
NFFT= 2^nextpow2(sum(inds2use_data));
[Ps_env, ~, bx1]= plot_dpss_psd(s_data_env(inds2use_data), fs_data, 'NW', nw, 'plotconf', true, 'nfft', NFFT);
[Psn_env, ~, bx2]= plot_dpss_psd(sn_data_env(inds2use_data), fs_data, 'NW', nw, 'plotconf', true, 'nfft', NFFT);
[Pnf_env, ~, bx3, Pnf_env_CI]= plot_dpss_psd([s_data_env(inds2use_nfloor) sn_data_env(inds2use_nfloor)], fs_data, 'NW', nw, 'plotconf', true, 'nfft', NFFT);
set(bx1(1), 'color', co(1,:), 'linestyle', '-', 'linew', lw3);
set(bx2(1), 'color', co(2,:), 'linestyle', '-', 'linew', lw3);
set(bx3(1), 'color', co(5,:), 'linestyle', '-', 'linew', lw2);
bx1(2).FaceColor= co(1,:);
bx2(2).FaceColor= co(2,:);
bx3(2).FaceColor= co(5,:);

Pxx= [Ps_env; Psn_env];
yl1= [max(Pxx)+5-yRange_ffr max(Pxx)+5];
% ylim(ylHard);
ylabel('');
ax_right(1)= gca;

lg= legend([ax bx1(1) bx2(1) bx3(1)], 'Stim', 'S-ENV', 'SN-ENV', 'NF', 'location', 'southwest');
set(lg, 'fontsize', 10);
set(gca, 'fontsize', fSize, 'xtick', xtick_vals, 'ytick', '', 'XTickLabel', xtick_labs);
title([ttlStr '-ENV']);
xlim(xlHard);

%%
subplot(nSProws, nSPcols, 2);
yyaxis left;
[Pxx, ~, ax]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
ylim([max(Pxx)+5-yRange_stim max(Pxx)+5]);
ylabel('');
ax_left(2)= gca;
set(gca, 'ytick', '');
set(ax, 'color', stimColor);
set(gca,'ycolor', stimColor);

yyaxis right;
hold on;
[Ps_tfs, ~, bx4]= plot_dpss_psd(s_data_tfs(inds2use_data) , fs_data, 'NW', nw, 'plotconf', true, 'nfft', NFFT);
[Psn_tfs, PSDfreq, bx5]= plot_dpss_psd(sn_data_tfs(inds2use_data) , fs_data, 'NW', nw, 'plotconf', true, 'nfft', NFFT);
[Pnf_tfs, ~, bx6, Pnf_tfs_CI]= plot_dpss_psd([s_data_tfs(inds2use_nfloor) sn_data_tfs(inds2use_nfloor)], fs_data, 'NW', nw, 'plotconf', true, 'nfft', NFFT);

set(bx4(1), 'color', co(1,:), 'linestyle', '-', 'linew', lw3);
set(bx5(1), 'color', co(2,:), 'linestyle', '-', 'linew', lw3);
set(bx6(1), 'color', co(5,:), 'linestyle', '-', 'linew', lw2);
bx4(2).FaceColor= co(1,:);
bx5(2).FaceColor= co(2,:);
bx6(2).FaceColor= co(5,:);


Pxx= [Ps_tfs; Psn_tfs];
yl2= [max(Pxx)+5-yRange_ffr max(Pxx)+5];
ax_right(2)= gca;
ylabel('FFR-PSD (dB/Hz)');

title([ttlStr '-TFS']);

lg= legend([ax bx4(1) bx5(1) bx6(1)], 'Stim', 'S-TFS', 'SN-TFS', 'NF', 'location', 'southwest');
set(lg, 'fontsize', 10);
set(gca, 'fontsize', fSize, 'xtick', xtick_vals, 'XTickLabel', xtick_labs);

%%
set(gcf, 'units', 'inches', 'position', [1 1 11 5]);
xlim(xlHard);

linkaxes(ax_left, 'y');
linkaxes(ax_right, 'y');
yl = [yl1; yl2];
% ylim([min(yl(:,1)) max(yl(:,2))]);
ylim(ylHard);

if saveFigs
    %     saveas(gcf, [fig_save_dir fName]);
    saveas(gcf, [fig_save_dir fName], 'png');
end

% Assign output
freqMin= 90;
freqMAX= 500;
valid_freq_inds= PSDfreq>freqMin & PSDfreq<freqMAX;

Ps_env_lin= db2mag(2*Ps_env); % factor of 2 because power to amplitude conversion
Psn_env_lin= db2mag(2*Psn_env);
Pnf_env_lin= db2mag(2*Pnf_env);
Pnf_env_CI_lin= db2mag(2*Pnf_env_CI);


Ps_tfs_lin= db2mag(2*Ps_tfs);
Psn_tfs_lin= db2mag(2*Psn_tfs);
Pnf_tfs_lin= db2mag(2*Pnf_tfs);
Pnf_tfs_CI_lin= db2mag(2*Pnf_tfs_CI);

PSD_struct.raw.ENV.S = .5*db(nansum(Ps_env_lin(valid_freq_inds)));
PSD_struct.raw.ENV.SN= .5*db(nansum(Psn_env_lin(valid_freq_inds)));
PSD_struct.raw.ENV.NF= .5*db(nansum(Pnf_env_lin(valid_freq_inds)));
PSD_struct.raw.ENV.NF_CI_low= .5*db(nansum(Pnf_env_CI_lin(valid_freq_inds,1)));
PSD_struct.raw.ENV.NF_CI_hi= .5*db(nansum(Pnf_env_CI_lin(valid_freq_inds,2)));

PSD_struct.raw.TFS.S = .5*db(nansum(Ps_tfs_lin(valid_freq_inds)));
PSD_struct.raw.TFS.SN= .5*db(nansum(Psn_tfs_lin(valid_freq_inds)));
PSD_struct.raw.TFS.NF= .5*db(nansum(Pnf_tfs_lin(valid_freq_inds)));
PSD_struct.raw.TFS.NF_CI_low= .5*db(nansum(Pnf_tfs_CI_lin(valid_freq_inds,1)));
PSD_struct.raw.TFS.NF_CI_hi= .5*db(nansum(Pnf_tfs_CI_lin(valid_freq_inds,2)));

PSD_struct.frac.ENV.S= db(nansum(Ps_env_lin(valid_freq_inds))/nansum(Ps_env_lin));
PSD_struct.frac.ENV.SN= db(nansum(Psn_env_lin(valid_freq_inds))/nansum(Psn_env_lin));
PSD_struct.frac.TFS.S= db(nansum(Ps_tfs_lin(valid_freq_inds))/nansum(Ps_tfs_lin));
PSD_struct.frac.TFS.SN= db(nansum(Psn_tfs_lin(valid_freq_inds))/nansum(Psn_tfs_lin));