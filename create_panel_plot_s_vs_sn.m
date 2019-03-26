function PSD_struct= create_panel_plot_s_vs_sn(figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, s_data_pos_filt, s_data_neg_filt, tStart, tEnd, fig_save_dir, fName, ttlStr, saveFigs)

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

yRange= 50;
lw= 2;
nw=2.5;
nSProws=1;
nSPcols=2;
fSize= 20;

figure(figHan);
clf;
co= get(gca, 'colororder');
lw3= 3;
% ylHard= [-50-yRange -50];
xlHard= [50 5e3];

%%
xticks= [10 100 1e3 4e3];
subplot(nSProws, nSPcols, 1);
yyaxis left;
[Pxx, ~, ~]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
ylim([max(Pxx)+5-yRange max(Pxx)+5]);
ax_left(1)= gca;
ylabel('Stim-PSD (dB/Hz)');

yyaxis right;
hold on;
[Ps_env, ~, bx1]= plot_dpss_psd(s_data_env(inds2use_data), fs_data, 'NW', nw);
[Psn_env, ~, bx2]= plot_dpss_psd(sn_data_env(inds2use_data), fs_data, 'NW', nw);
set(bx1, 'color', co(2,:), 'linestyle', '-', 'linew', lw3);
set(bx2, 'color', co(5,:), 'linestyle', '-', 'linew', lw3);
Pxx= [Ps_env; Psn_env];
yl1= [max(Pxx)+5-yRange max(Pxx)+5];
% ylim(ylHard);
ylabel('');
ax_right(1)= gca;

lg= legend('Stim', 'S-ENV', 'SN-ENV', 'location', 'southwest');
set(lg, 'fontsize', 10);
set(gca, 'fontsize', fSize, 'xtick', xticks, 'ytick', '');
title([ttlStr '-ENV']);
xlim(xlHard);

%%
subplot(nSProws, nSPcols, 2);
yyaxis left;
[Pxx, ~, ~]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
ylim([max(Pxx)+5-yRange max(Pxx)+5]);
ylabel('');
ax_left(2)= gca;
set(gca, 'ytick', '');

yyaxis right;
hold on;
[Ps_tfs, ~, bx3]= plot_dpss_psd(s_data_tfs(inds2use_data) , fs_data, 'NW', nw);
[Psn_tfs, PSDfreq, bx4]= plot_dpss_psd(sn_data_tfs(inds2use_data) , fs_data, 'NW', nw);
set(bx3, 'color', co(2,:), 'linestyle', '-', 'linew', lw3);
set(bx4, 'color', co(5,:), 'linestyle', '-', 'linew', lw3);
Pxx= [Ps_tfs; Psn_tfs];
yl2= [max(Pxx)+5-yRange max(Pxx)+5];
% ylim(ylHard);
ax_right(2)= gca;
ylabel('FFR-PSD (dB/Hz)');

title([ttlStr '-TFS']);

lg= legend('Stim', 'S-TFS', 'SN-TFS', 'location', 'southwest');
set(lg, 'fontsize', 10);
set(gca, 'fontsize', fSize, 'xtick', xticks);

%%
set(gcf, 'units', 'inches', 'position', [1 1 10 4]);
xlim(xlHard);

linkaxes(ax_left, 'y');
linkaxes(ax_right, 'y');
yl = [yl1; yl2];
ylim([min(yl(:,1)) max(yl(:,2))]);

if saveFigs
%     saveas(gcf, [fig_save_dir fName]);
    saveas(gcf, [fig_save_dir fName], 'png');
end

% Assign output
freqMAX= 500;
valid_freq_inds= PSDfreq<freqMAX;

Ps_env_lin= db2mag(2*Ps_env); % factor of 2 because power to amplitude conversion
Psn_env_lin= db2mag(2*Psn_env);
Ps_tfs_lin= db2mag(2*Ps_tfs);
Psn_tfs_lin= db2mag(2*Psn_tfs);

PSD_struct.raw.ENV.S= .5*db(sum(Ps_env_lin(valid_freq_inds)));
PSD_struct.raw.ENV.SN= .5*db(sum(Psn_env_lin(valid_freq_inds)));
PSD_struct.raw.TFS.S= .5*db(sum(Ps_tfs_lin(valid_freq_inds)));
PSD_struct.raw.TFS.SN= .5*db(sum(Psn_tfs_lin(valid_freq_inds)));

PSD_struct.frac.ENV.S= db(sum(Ps_env_lin(valid_freq_inds))/sum(Ps_env_lin));
PSD_struct.frac.ENV.SN= db(sum(Psn_env_lin(valid_freq_inds))/sum(Psn_env_lin));
PSD_struct.frac.TFS.S= db(sum(Ps_tfs_lin(valid_freq_inds))/sum(Ps_tfs_lin));
PSD_struct.frac.TFS.SN= db(sum(Psn_tfs_lin(valid_freq_inds))/sum(Psn_tfs_lin));