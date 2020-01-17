function PSD_struct= create_panel_plot_s_vs_sn_tdt...
    (figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, s_data_pos_filt, s_data_neg_filt, ...
    sn_nf_pos_filt, sn_nf_neg_filt, s_nf_pos_filt, s_nf_neg_filt, ...
    tStart, tEnd, fig_save_dir, fName, plotNF, saveFigs)

sn_data_env= (sn_data_pos_filt+sn_data_neg_filt)/2;
sn_data_tfs= (sn_data_pos_filt-sn_data_neg_filt)/2;

s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;

sn_nf_env= (sn_nf_pos_filt+sn_nf_neg_filt)/2;
sn_nf_tfs= (sn_nf_pos_filt-sn_nf_neg_filt)/2;

s_nf_env= (s_nf_pos_filt+s_nf_neg_filt)/2;
s_nf_tfs= (s_nf_pos_filt-s_nf_neg_filt)/2;


if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end

t_sig= (1:length(sig))/fs_sig;
t_data= (1:length(sn_data_pos_filt))/fs_data;

inds2use_stim= t_sig>tStart & t_sig<tEnd;
inds2use_data= t_data>tStart & t_data<tEnd;

yRange_stim= 50;
lw2= 2;
nw=7;
nw_nf=nw;
nSProws=1;
nSPcols=2;
fSize= 20;

figure(figHan);
clf;

lw3= 3;
ylHard= [-50-yRange_stim -50];
xlHard= [75 750];
stimColor= .4*[1 1 1];

nfft= 2^nextpow2(sum(inds2use_data));

%%
xtick_vals= [10 100 500 750];
tick_len= [.025 .025];
xtick_labs= cellfun(@(x) num2str(x), num2cell(xtick_vals), 'UniformOutput', false);
subplot(nSProws, nSPcols, 1);
co= set_colblind_order();

yyaxis left;
[Pxx, ~, ax]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
ylim([max(Pxx)+5-yRange_stim max(Pxx)+5]);
ax_left(1)= gca;
ylabel('Stim-PSD (dB/Hz)');
set(ax, 'color', stimColor);
set(gca,'ycolor', 'k', 'box', 'off', 'LineWidth', 1.5, 'TickLength', tick_len);
grid off;

yyaxis right;
hold on;
[Ps_env, ~, bx1]= plot_dpss_psd(s_data_env(inds2use_data), fs_data, 'NW', nw, 'plotconf', true, 'nfft', nfft);
[Psn_env, ~, bx2]= plot_dpss_psd(sn_data_env(inds2use_data), fs_data, 'NW', nw, 'plotconf', true, 'nfft', nfft);
grid off;

if plotNF
    [Ps_nf_env, ~, ~, ~]= plot_dpss_psd(s_nf_env(inds2use_data), fs_data, 'NW', nw_nf, 'plot', false, 'nfft', nfft);
    [Psn_nf_env, freq_nf, ~, Pnf_env_CI]= plot_dpss_psd(sn_nf_env(inds2use_data), fs_data, 'NW', nw_nf, 'plot', false, 'nfft', nfft);
    bx3(1)= plot(freq_nf, Psn_nf_env, '-');
    cx= fill([freq_nf;flipud(freq_nf)], [nanmean(Pnf_env_CI(:, 1:2:end) ,2); flipud(nanmean(Pnf_env_CI(:, 2:2:end) ,2))], get(bx3(1), 'color'));
    bx3= [bx3; cx];
    
    set(bx3(1), 'color', co(5,:), 'linestyle', '-', 'linew', lw2);
    bx3(2).FaceColor= co(5,:);
    bx3(2).FaceAlpha= .4;
    bx3(2).EdgeAlpha= 0;
    
    lg= legend([ax bx1(1) bx2(1) bx3(1)], 'Stim', 'S-ENV', 'SN-ENV', 'NF', 'location', 'southwest', 'box', 'off');
else
    lg= legend([ax bx1(1) bx2(1) ], 'Stim', 'S-ENV', 'SN-ENV', 'location', 'southwest', 'box', 'off');
end

set(bx1(1), 'color', co(1,:), 'linestyle', '-', 'linew', lw3);
set(bx2(1), 'color', co(2,:), 'linestyle', '-', 'linew', lw3);
bx1(2).FaceColor= co(1,:);
bx2(2).FaceColor= co(2,:);

ylabel('');
ax_right(1)= gca;

set(lg, 'fontsize', fSize);
set(gca, 'fontsize', fSize, 'xtick', xtick_vals, 'ytick', '', 'XTickLabel', xtick_labs);
% title([ttlStr '-ENV']);
title('ENV')
xlim(xlHard);

%%
subplot(nSProws, nSPcols, 2);
set_colblind_order();
yyaxis left;
[Pxx, ~, ax]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
ylim([max(Pxx)+5-yRange_stim max(Pxx)+5]);
ylabel('');
ax_left(2)= gca;
set(gca, 'ytick', '');
set(ax, 'color', stimColor);
set(gca,'ycolor', stimColor, 'box', 'off');
grid off;

yyaxis right;
hold on;
[Ps_tfs, ~, bx4]= plot_dpss_psd(s_data_tfs(inds2use_data) , fs_data, 'NW', nw, 'plotconf', true, 'nfft', nfft);
[Psn_tfs, freq_data, bx5]= plot_dpss_psd(sn_data_tfs(inds2use_data) , fs_data, 'NW', nw, 'plotconf', true, 'nfft', nfft);

if plotNF
    [Ps_nf_tfs, ~, ~, ~]= plot_dpss_psd(s_nf_tfs(inds2use_data), fs_data, 'NW', nw_nf, 'plot', false, 'nfft', nfft);
    [Psn_nf_tfs, freq_nf, ~, Pnf_tfs_CI]= plot_dpss_psd(sn_nf_tfs(inds2use_data), fs_data, 'NW', nw_nf, 'plot', false, 'nfft', nfft);
    bx6(1)= plot(freq_nf, Psn_nf_tfs, '-');
    cx= fill([freq_nf;flipud(freq_nf)], [nanmean(Pnf_tfs_CI(:, 1:2:end) ,2); flipud(nanmean(Pnf_tfs_CI(:, 2:2:end) ,2))], get(bx6(1), 'color'));
    bx6= [bx6; cx];
    set(bx6(1), 'color', co(5,:), 'linestyle', '-', 'linew', lw2);
    bx6(2).FaceColor= co(5,:);
    bx6(2).FaceAlpha= .4;
    bx6(2).EdgeAlpha= 0;
    lg= legend([ax bx4(1) bx5(1) bx6(1)], 'Stim', 'S-TFS', 'SN-TFS', 'NF', 'location', 'southwest', 'box', 'off');
else
    lg= legend([ax bx4(1) bx5(1)], 'Stim', 'S-TFS', 'SN-TFS', 'location', 'southwest', 'box', 'off');
    freq_nf= freq_data;
    Ps_nf_env= 0*Ps_tfs;
    Psn_nf_env= 0*Ps_tfs;
    Pnf_env_CI= 0*[Ps_tfs(:), Ps_tfs(:)];
    Ps_nf_tfs= 0*Ps_tfs;
    Psn_nf_tfs= 0*Ps_tfs;
    Pnf_tfs_CI= 0*[Ps_tfs(:), Ps_tfs(:)];
end

set(bx4(1), 'color', co(1,:), 'linestyle', '-', 'linew', lw3);
set(bx5(1), 'color', co(2,:), 'linestyle', '-', 'linew', lw3);
bx4(2).FaceColor= co(1,:);
bx5(2).FaceColor= co(2,:);

ax_right(2)= gca;
ylabel('FFR-PSD (dB/Hz)');

% title([ttlStr '-TFS']);
title('TFS')

set(lg, 'fontsize', fSize);
set(gca, 'fontsize', fSize, 'xtick', xtick_vals, 'XTickLabel', xtick_labs, 'LineWidth', 1.5, 'TickLength', tick_len);
grid off;

%%
set(gcf, 'units', 'inches', 'position', [1 1 11 5],'Renderer','painters');

linkaxes(ax_left, 'y');
linkaxes(ax_right, 'y');
% yl = [yl1; yl2];
% ylim([min(yl(:,1)) max(yl(:,2))]);
ylim(ylHard);
xlim(xlHard);

add_subplot_letter(1, 2, 30);
if saveFigs
    %     saveas(gcf, [fig_save_dir fName]);
    saveas(gcf, [fig_save_dir fName], 'png');
    saveas(gcf, [fig_save_dir fName], 'tiff');
end

%% Assign output
freqMin= 90;
freqMAX= 500;

valid_freq_inds_data= freq_data>freqMin & freq_data<freqMAX;
valid_freq_inds_nf= freq_nf>freqMin & freq_nf<freqMAX;

Ps_env_lin= db2mag(2*Ps_env); % factor of 2 because power to amplitude conversion
Psn_env_lin= db2mag(2*Psn_env);
Ps_nf_env_lin= db2mag(2*Ps_nf_env);
Psn_nf_env_lin= db2mag(2*Psn_nf_env);
Psn_nf_env_CI_lin= db2mag(2*Pnf_env_CI);


Ps_tfs_lin= db2mag(2*Ps_tfs);
Psn_tfs_lin= db2mag(2*Psn_tfs);
Ps_nf_tfs_lin= db2mag(2*Ps_nf_tfs);
Psn_nf_tfs_lin= db2mag(2*Psn_nf_tfs);
Pnf_tfs_CI_lin= db2mag(2*Pnf_tfs_CI);

data_scale_factor= fs_data/nfft;
nf_scale_factor= fs_data/nfft;


PSD_struct.raw.ENV.S = .5*db(nansum(Ps_env_lin(valid_freq_inds_data))*data_scale_factor);
PSD_struct.raw.ENV.SN= .5*db(nansum(Psn_env_lin(valid_freq_inds_data))*data_scale_factor);
PSD_struct.raw.ENV.NF_S= .5*db(nansum(Ps_nf_env_lin(valid_freq_inds_nf))*nf_scale_factor);
PSD_struct.raw.ENV.NF_SN= .5*db(nansum(Psn_nf_env_lin(valid_freq_inds_nf))*nf_scale_factor);
PSD_struct.raw.ENV.NF_SN_CI_low= .5*db(nansum(Psn_nf_env_CI_lin(valid_freq_inds_nf,1))*nf_scale_factor);
PSD_struct.raw.ENV.NF_SN_CI_hi= .5*db(nansum(Psn_nf_env_CI_lin(valid_freq_inds_nf,2))*nf_scale_factor);

PSD_struct.raw.TFS.S = .5*db(nansum(Ps_tfs_lin(valid_freq_inds_data))*data_scale_factor);
PSD_struct.raw.TFS.SN= .5*db(nansum(Psn_tfs_lin(valid_freq_inds_data))*data_scale_factor);
PSD_struct.raw.TFS.NF_S= .5*db(nansum(Ps_nf_tfs_lin(valid_freq_inds_nf))*nf_scale_factor);
PSD_struct.raw.TFS.NF_SN= .5*db(nansum(Psn_nf_tfs_lin(valid_freq_inds_nf))*nf_scale_factor);
PSD_struct.raw.TFS.NF_SN_CI_low= .5*db(nansum(Pnf_tfs_CI_lin(valid_freq_inds_nf,1))*nf_scale_factor);
PSD_struct.raw.TFS.NF_SN_CI_hi= .5*db(nansum(Pnf_tfs_CI_lin(valid_freq_inds_nf,2))*nf_scale_factor);

PSD_struct.frac.ENV.S= db(nansum(Ps_env_lin(valid_freq_inds_data))/nansum(Ps_env_lin));
PSD_struct.frac.ENV.SN= db(nansum(Psn_env_lin(valid_freq_inds_data))/nansum(Psn_env_lin));
PSD_struct.frac.TFS.S= db(nansum(Ps_tfs_lin(valid_freq_inds_data))/nansum(Ps_tfs_lin));
PSD_struct.frac.TFS.SN= db(nansum(Psn_tfs_lin(valid_freq_inds_data))/nansum(Psn_tfs_lin));