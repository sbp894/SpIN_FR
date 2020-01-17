function legHan= create_panel_plot_s_vs_sn_tdt_at2snrs...
    (figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, s_data_pos_filt, s_data_neg_filt, ...
    tStart, tEnd, plot_S_sig)

sn_data_env= (sn_data_pos_filt+sn_data_neg_filt)/2;
sn_data_tfs= (sn_data_pos_filt-sn_data_neg_filt)/2;

s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;

t_sig= (1:length(sig))/fs_sig;
t_data= (1:length(sn_data_pos_filt))/fs_data;

inds2use_stim= t_sig>tStart & t_sig<tEnd;
inds2use_data= t_data>tStart & t_data<tEnd;

yRange_stim= 50;
nw=7;
nSProws=1;
nSPcols=2;
fSize= 20;

figure(figHan);
% clf;

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
% co= set_colblind_order();
% co= co([1 2 6], :);
co= [get_color('b'); get_color('r'); get_color('g')];

yyaxis left;
if plot_S_sig
    [Pxx, ~, ax]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
    ylim([max(Pxx)+5-yRange_stim max(Pxx)+5]);
    ax_left(1)= gca;
    ylabel('Stim-PSD (dB/Hz)');
    set(ax, 'color', stimColor);
    set(gca,'ycolor', 'k', 'box', 'off', 'LineWidth', 1.5, 'TickLength', tick_len);
    grid off;
    legHan= nan(3,1);
    legHan(1)= ax;
else 
    legHan= nan;
end

yyaxis right;
hold on;
if plot_S_sig
    [~, ~, bx1]= plot_dpss_psd(s_data_env(inds2use_data), fs_data, 'NW', nw, 'plotconf', true, 'nfft', nfft);
    set(bx1(1), 'color', co(1,:), 'linestyle', '-', 'linew', lw3);
    bx1(2).FaceColor= co(1,:);
    legHan(2)= bx1(1);
end
[~, ~, bx2]= plot_dpss_psd(sn_data_env(inds2use_data), fs_data, 'NW', nw, 'plotconf', true, 'nfft', nfft);
set(bx2(1), 'color', co(2+plot_S_sig,:), 'linestyle', '-', 'linew', lw3);
bx2(2).FaceColor= co(2+plot_S_sig,:);
grid off;
legHan(end)= bx2(1);

% lg= legend([ax bx1(1) bx2(1) ], 'Stim', 'S-ENV', 'SN-ENV', 'location', 'southwest', 'box', 'off');


ylabel('');
ax_right(1)= gca;

% set(lg, 'fontsize', fSize);
set(gca, 'fontsize', fSize, 'xtick', xtick_vals, 'ytick', '', 'XTickLabel', xtick_labs);
% title([ttlStr '-ENV']);
title('ENV')
xlim(xlHard);

%%
subplot(nSProws, nSPcols, 2);
set_colblind_order();
if plot_S_sig
    yyaxis left;
    [Pxx, ~, ax]= plot_dpss_psd(sig(inds2use_stim), fs_sig, 'NW', nw);
    ylim([max(Pxx)+5-yRange_stim max(Pxx)+5]);
    ylabel('');
    ax_left(2)= gca;
    set(gca, 'ytick', '');
    set(ax, 'color', stimColor);
    set(gca,'ycolor', stimColor, 'box', 'off');
    grid off;
    linkaxes(ax_left, 'y');
end

yyaxis right;
hold on;
if plot_S_sig
    [~, ~, bx4]= plot_dpss_psd(s_data_tfs(inds2use_data) , fs_data, 'NW', nw, 'plotconf', true, 'nfft', nfft);
    set(bx4(1), 'color', co(1,:), 'linestyle', '-', 'linew', lw3);
    bx4(2).FaceColor= co(1,:);
end
[~, ~, bx5]= plot_dpss_psd(sn_data_tfs(inds2use_data) , fs_data, 'NW', nw, 'plotconf', true, 'nfft', nfft);
set(bx5(1), 'color', co(2+plot_S_sig,:), 'linestyle', '-', 'linew', lw3);
bx5(2).FaceColor= co(2+plot_S_sig,:);
xlabel('');


% lg= legend([ax bx4(1) bx5(1)], 'Stim', 'S-TFS', 'SN-TFS', 'location', 'southwest', 'box', 'off');



ax_right(2)= gca;
linkaxes(ax_right, 'y');

ylabel('FFR-PSD (dB/Hz)');

% title([ttlStr '-TFS']);
title('TFS')

% set(lg, 'fontsize', fSize);
set(gca, 'fontsize', fSize, 'xtick', xtick_vals, 'XTickLabel', xtick_labs, 'LineWidth', 1.5, 'TickLength', tick_len);
grid off;

ylim(ylHard);
xlim(xlHard);

add_subplot_letter(1, 2, 30);
