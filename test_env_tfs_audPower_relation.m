function test_env_tfs_audPower_relation...
    (pool_lf_power_audio, pool_hf_power_audio, pool_env_power_ffr, pool_tfs_power_ffr, nhInds, hiInds, LatexDir, postFix, data_save_dir, nh_ChinIDs, hi_ChinIDs)

ValidRows= any(~isnan(pool_lf_power_audio), 2);
notNANinds= ~isnan(pool_lf_power_audio(find(ValidRows, 1, 'last' ), :));

if ~all(ValidRows)
    error('NH and HI inds may be shifted ');
end

tick_len= [.025 .025];
pool_lf_power_audio_db= dbspl(sqrt(pool_lf_power_audio(ValidRows, notNANinds)));
pool_hf_power_audio_db= dbspl(sqrt(pool_hf_power_audio(ValidRows, notNANinds)));
pool_env_power_ffr_db= db(sqrt(pool_env_power_ffr(ValidRows, notNANinds)));
pool_tfs_power_ffr_db= db(sqrt(pool_tfs_power_ffr(ValidRows, notNANinds)));

set_maxFFR_val_to_zero= 1;
if set_maxFFR_val_to_zero
    pool_env_power_ffr_db= pool_env_power_ffr_db-max(pool_env_power_ffr_db(:));
    pool_tfs_power_ffr_db= pool_tfs_power_ffr_db-max(pool_tfs_power_ffr_db(:));
end

ignore_high_LFpower= 0;
if ignore_high_LFpower
    power_Thresh= 66;
    pool_lf_power_audio_db(pool_lf_power_audio_db>power_Thresh)= nan;
end


pValThresh= 1e-3;

params_nh= struct('x_txt_val', .92, 'y_txt_val', .45, 'y_txt_gap', .08, 'fSize', 16, 'pValThresh', pValThresh, 'title', 'NH');
params_hi= struct('x_txt_val', .92, 'y_txt_val', .95, 'y_txt_gap', .08, 'fSize', 16, 'pValThresh', pValThresh, 'title', 'HI');


mrkSize= 7;
lw2= 1.5;
lw3= 5;
fSize= 20;


xtick_lf= 35:5:75;
xtick_hf= 35:10:75;

figure(88);
clf;
co= set_colblind_order();

for chinVar = 1:size(pool_lf_power_audio_db,1)
    enterFlag= 0;
    if ismember(chinVar, nhInds)
        enterFlag= 1;
        plotColor= co(1,:);
        markerType= 's';
    elseif ismember(chinVar, hiInds)
        enterFlag= 1;
        plotColor= co(2,:);
        markerType= 'd';
    end
    
    if enterFlag
        % -----------------
        ax(1)= subplot(221);
        hold on;
        plot(pool_lf_power_audio_db(chinVar,:), pool_env_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', mrkSize, 'LineWidth', lw2);
        ylabel('$ENV_{FFR}^{power}$ (dB)', 'interpreter', 'latex');
        
        
        % -----------------
        ax(2)= subplot(222);
        hold on;
        plot(pool_hf_power_audio_db(chinVar,:), pool_env_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', mrkSize, 'LineWidth', lw2);
        
        % -----------------
        ax(3)= subplot(223);
        hold on;
        plot(pool_lf_power_audio_db(chinVar,:), pool_tfs_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', mrkSize, 'LineWidth', lw2);
        xlabel('$LF_{stimulus}^{power}$ (dB SPL)', 'interpreter', 'latex');
        ylabel('$TFS_{FFR}^{power}$ (dB)', 'interpreter', 'latex');
        
        % -----------------
        ax(4)= subplot(224);
        hold on;
        plot(pool_hf_power_audio_db(chinVar,:), pool_tfs_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', mrkSize, 'LineWidth', lw2);
        xlabel('$HF_{stimulus}^{power}$ (dB SPL)', 'interpreter', 'latex');
    end
end

nh_lf_aud= pool_lf_power_audio_db(nhInds,:); % values for LF are (and should be) same for all animals
nh_lf_aud= nh_lf_aud(:);
nh_hf_aud= pool_hf_power_audio_db(nhInds,:); % similarly, values for HF are (and should be) same for all animals
nh_hf_aud= nh_hf_aud(:);
hi_lf_aud= pool_lf_power_audio_db(hiInds,:);
hi_lf_aud= hi_lf_aud(:);
hi_hf_aud= pool_hf_power_audio_db(hiInds,:);
hi_hf_aud= hi_hf_aud(:);

nh_env_ffr= pool_env_power_ffr_db(nhInds,:);
nh_env_ffr= nh_env_ffr(:);
nh_tfs_ffr= pool_tfs_power_ffr_db(nhInds,:);
nh_tfs_ffr= nh_tfs_ffr(:);
hi_env_ffr= pool_env_power_ffr_db(hiInds,:);
hi_env_ffr= hi_env_ffr(:);
hi_tfs_ffr= pool_tfs_power_ffr_db(hiInds,:);
hi_tfs_ffr= hi_tfs_ffr(:);


lf_x_vals= pool_lf_power_audio_db(nhInds(1),:)';
hf_x_vals= pool_hf_power_audio_db(nhInds(1),:)';

subplot(221);
mdl1_nh= fitlm(nh_lf_aud(:), nh_env_ffr(:));
nh_env_lf_fit= predict(mdl1_nh, lf_x_vals);
mdl1_hi= fitlm(hi_lf_aud(:), hi_env_ffr(:));
hi_env_lf_fit= predict(mdl1_hi, lf_x_vals);
plot(lf_x_vals, nh_env_lf_fit, 'color', co(1,:), 'LineWidth', lw3);
plot(lf_x_vals, hi_env_lf_fit, 'color', co(2,:), 'LineWidth', lw3);
set(gca, 'FontSize', fSize, 'XTick', xtick_lf, 'LineWidth', 1.5, 'TickLength', tick_len);
grid off;

add_stat_txt(mdl1_hi, params_hi);
add_stat_txt(mdl1_nh, params_nh);

subplot(222);
mdl2_nh= fitlm(nh_hf_aud, nh_env_ffr);
nh_env_hf_fit= predict(mdl2_nh, hf_x_vals);
mdl2_hi= fitlm(hi_hf_aud, hi_env_ffr);
hi_env_hf_fit= predict(mdl2_hi, hf_x_vals);
plot(hf_x_vals, nh_env_hf_fit, 'color', co(1,:), 'LineWidth', lw3);
plot(hf_x_vals, hi_env_hf_fit, 'color', co(2,:), 'LineWidth', lw3);
set(gca, 'FontSize', fSize, 'XTick', xtick_hf, 'LineWidth', 1.5, 'TickLength', tick_len);
grid off;

add_stat_txt(mdl2_hi, params_hi);
add_stat_txt(mdl2_nh, params_nh);


subplot(223);
mdl3_nh= fitlm(nh_lf_aud(:), nh_tfs_ffr(:));
nh_tfs_lf_fit= predict(mdl3_nh, lf_x_vals);
mdl3_hi= fitlm(hi_lf_aud(:), hi_tfs_ffr(:));
hi_tfs_lf_fit= predict(mdl3_hi, lf_x_vals);
plot(lf_x_vals, nh_tfs_lf_fit, 'color', co(1,:), 'LineWidth', lw3);
plot(lf_x_vals, hi_tfs_lf_fit, 'color', co(2,:), 'LineWidth', lw3);
set(gca, 'FontSize', fSize, 'XTick', xtick_lf, 'LineWidth', 1.5, 'TickLength', tick_len);
grid off;

add_stat_txt(mdl3_hi, params_hi);
add_stat_txt(mdl3_nh, params_nh);


subplot(224);
mdl4_nh= fitlm(nh_hf_aud, nh_tfs_ffr);
nh_env_hf_fit= predict(mdl4_nh, hf_x_vals);
mdl4_hi= fitlm(hi_hf_aud, hi_tfs_ffr);
hi_env_hf_fit= predict(mdl4_hi, hf_x_vals);
plot(hf_x_vals, nh_env_hf_fit, 'color', co(1,:), 'LineWidth', lw3);
plot(hf_x_vals, hi_env_hf_fit, 'color', co(2,:), 'LineWidth', lw3);
set(gca, 'FontSize', fSize, 'XTick', xtick_hf, 'LineWidth', 1.5, 'TickLength', tick_len);
grid off;

add_stat_txt(mdl4_hi, params_hi);
add_stat_txt(mdl4_nh, params_nh);


% mdl_lf_hf= fitlm(lf_x_vals, hf_x_vals)

linkaxes(ax([1 3]), 'x')
subplot(221);
xlim([min(pool_lf_power_audio_db(:))-1 max(pool_lf_power_audio_db(:))+1.5]);
linkaxes(ax([2 4]), 'x')
subplot(222);
xlim([min(pool_hf_power_audio_db(:))-1 max(pool_hf_power_audio_db(:))+4]);

linkaxes(ax([1 2]), 'y')
ylim([-18 0]);
linkaxes(ax([3 4]), 'y')


subplot(221);
ll(1)= plot(nan, nan, '-s', 'color', co(1,:), 'MarkerSize', mrkSize, 'LineWidth', lw2);
ll(2)= plot(nan, nan, '-d', 'color', co(2,:), 'MarkerSize', mrkSize, 'LineWidth', lw2);
lg= legend(ll, 'NH', 'HI', 'Location', 'south', 'box', 'off', 'Orientation', 'horizontal');
lg.FontSize= fSize;


set(gcf, 'Units', 'inches', 'Position', [1 1 11 8]);


add_subplot_letter(2, 2, 30, 0, 1.09);

fName= ['Figure_Out/all_corr_sfr_4panel' postFix];
if ignore_high_LFpower
    fName= [fName '_ignoredHF'];
end
saveas(gcf, fName, 'png');

fNameLatex= [LatexDir 'all_corr_sfr_4panel' postFix];
saveas(gcf, fNameLatex, 'epsc');

save([data_save_dir 'nh_hi_data_lf_tfs.mat'], 'pool_lf_power_audio_db', 'pool_hf_power_audio_db', 'pool_tfs_power_ffr_db', 'pool_env_power_ffr_db', ...
    'nhInds', 'hiInds', 'nh_ChinIDs', 'hi_ChinIDs');