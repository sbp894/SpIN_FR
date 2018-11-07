function create_panel_plot_include_hilbert(figHan, fs_sig, sig, fs_data, data_pos_filt, data_neg_filt, data_env, data_tfs, tStart, tEnd, fig_save_dir, fName, ttlStr, saveFigs)

t_sig= (1:length(sig))/fs_sig;
t_data= (1:length(data_pos_filt))/fs_data;

inds2use_stim= t_sig>tStart & t_sig<tEnd;
inds2use_data= t_data>tStart & t_data<tEnd;

Ashift= 2*max(abs([data_pos_filt(inds2use_data) data_neg_filt(inds2use_data)]));

yRange= 100;
lw= 2;
nw=2.5;
nSProws=2;
nSPcols=2;

figure(figHan);
clf;

%%
subplot(nSProws, nSPcols, 1:2);

fSize= 16;
ax(1)= plot(t_data, 2*Ashift + data_pos_filt, 'linew', lw);
hold on;
ax(2)= plot(t_data, 2*Ashift + data_neg_filt, 'linew', lw);
ax(3)= plot(t_data, Ashift + abs(hilbert(data_pos_filt)), 'linew', lw);
ax(4)= plot(t_data, Ashift + abs(hilbert(data_neg_filt)), 'linew', lw);
ax(5)= plot(t_data, data_env, 'linew', lw);
ax(6)= plot(t_data,- Ashift + data_tfs, 'linew', lw);
ax(7)= plot(t_sig, - 2*Ashift + sig*Ashift, 'k', 'linew', lw);
grid on;
lg= legend(ax, '+ve pol', '-ve pol', '+ve Hilb', '-ve Hilb', 'sum of pols', 'diff of pols',  'stim');
set(lg, 'fontsize', 10, 'location', 'eastoutside');
set(gca, 'fontsize', fSize);
xlabel('time (sec)');
ylabel('Amp');
title(ttlStr);
xlim([tStart-5e-3 tEnd+5e-3]);

%%
subplot(nSProws, nSPcols, 3);
yyaxis left;
[Pxx, ~, ~]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
ylim([max(Pxx)+5-yRange max(Pxx)+5]);

yyaxis right;
hold on;
[Pxx_env, ~, bx1]= plot_dpss_psd(data_env(inds2use_data), fs_data, 'NW', nw);
[Pxx_tfs, ~, bx2]= plot_dpss_psd(data_tfs(inds2use_data), fs_data, 'NW', nw);
set(bx1, 'color', get(ax(5), 'color'), 'linestyle', '-');
set(bx2, 'color', get(ax(6), 'color'), 'linestyle', '-');
Pxx= [Pxx_env; Pxx_tfs];
ylim([max(Pxx)+5-yRange max(Pxx)+5]);

legend('Stim', 'Data-ENV', 'Data-TFS');
set(gca, 'fontsize', fSize);

%%
subplot(nSProws, nSPcols, 4);
yyaxis left;
[Pxx, ~, ~]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
ylim([max(Pxx)+5-yRange max(Pxx)+5]);

yyaxis right;
hold on;
[Pxx_env, ~, bx3]= plot_dpss_psd(abs(hilbert(data_pos_filt(inds2use_data))), fs_data, 'NW', nw);
[Pxx_tfs, ~, bx4]= plot_dpss_psd(abs(hilbert(data_neg_filt(inds2use_data))), fs_data, 'NW', nw);
set(bx3, 'color', get(ax(3), 'color'), 'linestyle', '-');
set(bx4, 'color', get(ax(4), 'color'), 'linestyle', '-');
Pxx= [Pxx_env; Pxx_tfs];
ylim([max(Pxx)+5-yRange max(Pxx)+5]);

legend('Stim', 'Hilb (+ve)', 'Hilb (-ve)');
set(gca, 'fontsize', fSize);

%%
set(gcf, 'units', 'normalized', 'position', [0 0 1 1]);

if saveFigs
    saveas(gcf, [fig_save_dir fName]);
    saveas(gcf, [fig_save_dir fName], 'tiff');
end