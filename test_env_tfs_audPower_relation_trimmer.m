function test_env_tfs_audPower_relation_trimmer(pool_lf_power_audio, pool_hf_power_audio, pool_env_power_ffr, pool_tfs_power_ffr, nhInds, hiInds, LatexDir)

pool_lf_power_audio_db= dbspl(sqrt(pool_lf_power_audio));
pool_hf_power_audio_db= dbspl(sqrt(pool_hf_power_audio));
pool_env_power_ffr_db= db(sqrt(pool_env_power_ffr));
pool_tfs_power_ffr_db= db(sqrt(pool_tfs_power_ffr));

set_max_val_to_zero= 1;
if set_max_val_to_zero
    pool_env_power_ffr_db= pool_env_power_ffr_db-max(pool_env_power_ffr_db(:));
    pool_tfs_power_ffr_db= pool_tfs_power_ffr_db-max(pool_tfs_power_ffr_db(:));
end

mrkSize= 8;
lw2= 1.5;
lw3= 5;
fSize= 18;
fSize_txt= 16;
xtick_lf= 35:5:75;
xtick_hf= 35:10:75;


figure(88);
clf;
co= get(gca, 'ColorOrder');
for chinVar = 1:size(pool_lf_power_audio_db,1)
    enterFlag= 0;
    if ismember(chinVar, nhInds)
        enterFlag= 1;
        plotColor= co(1,:);
        markerType= 'o';
    elseif ismember(chinVar, hiInds)
        enterFlag= 1;
        plotColor= co(2,:);
        markerType= 'd';
    end
    
    if enterFlag
        
        % -----------------
        ax(1)= subplot(211);
        hold on;
        plot(pool_hf_power_audio_db(chinVar,:), pool_env_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', mrkSize, 'LineWidth', lw2);
        xlabel('$HF_{stimulus} (dB)$', 'Interpreter', 'latex');
        ylabel('$ENV_{FFR} (dB)$', 'Interpreter', 'latex');
        
        % -----------------
        ax(2)= subplot(212);
        hold on;
        plot(pool_lf_power_audio_db(chinVar,:), pool_tfs_power_ffr_db(chinVar,:), markerType, 'color', plotColor, 'MarkerSize', mrkSize, 'LineWidth', lw2);
        xlabel('$LF_{stimulus}$ (dB)', 'Interpreter', 'latex');
        ylabel('$TFS_{FFR}$ (dB)', 'Interpreter', 'latex');
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


subplot(211);
mdl2_nh= fitlm(nh_hf_aud, nh_env_ffr)
nh_env_hf_fit= predict(mdl2_nh, hf_x_vals);
mdl2_hi= fitlm(hi_hf_aud, hi_env_ffr)
hi_env_hf_fit= predict(mdl2_hi, hf_x_vals);
plot(hf_x_vals, nh_env_hf_fit, 'b', 'LineWidth', lw3);
plot(hf_x_vals, hi_env_hf_fit, 'r', 'LineWidth', lw3);
set(gca, 'FontSize', fSize, 'XTick', xtick_hf);
xlim([min(pool_hf_power_audio_db(:))-1 max(pool_hf_power_audio_db(:))+1]);
grid on;

pValThresh= 1e-4;
if mdl2_nh.Coefficients.pValue(2)>pValThresh
    text(.1,.1,sprintf('$p=%.2f, R^2=%.2f$', mdl2_nh.Coefficients.pValue(2), mdl2_nh.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold', 'Color', co(1,:));
else
    text(.1,.1,sprintf('$p<%.4f, R^2=%.2f$', pValThresh, mdl2_nh.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold', 'Color', co(1,:));
end
if mdl2_hi.Coefficients.pValue(2)>pValThresh
    text(.55,.1,sprintf('$p=%.2f, R^2=%.2f$', mdl2_hi.Coefficients.pValue(2), mdl2_hi.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold', 'Color', co(2,:));
else
    text(.55,.1,sprintf('$p<%.4f, R^2=%.2f$', pValThresh, mdl2_hi.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold', 'Color', co(2,:));
end

subplot(212);
mdl3_nh= fitlm(nh_lf_aud(:), nh_tfs_ffr(:))
nh_tfs_lf_fit= predict(mdl3_nh, lf_x_vals);
mdl3_hi= fitlm(hi_lf_aud(:), hi_tfs_ffr(:))
hi_tfs_lf_fit= predict(mdl3_hi, lf_x_vals);
plot(lf_x_vals, nh_tfs_lf_fit, 'b', 'LineWidth', lw3);
plot(lf_x_vals, hi_tfs_lf_fit, 'r', 'LineWidth', lw3);
set(gca, 'FontSize', fSize, 'XTick', xtick_lf);
grid on;
xlim([min(pool_lf_power_audio_db(:))-1 max(pool_lf_power_audio_db(:))+1]);

if mdl3_nh.Coefficients.pValue(2)>pValThresh
    text(.1,.05,sprintf('$p=%.2f, R^2=%.2f$', mdl3_nh.Coefficients.pValue(2), mdl3_nh.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold', 'Color', co(1,:));
else
    text(.1,.05,sprintf('$p<%.4f, R^2=%.2f$', pValThresh, mdl3_nh.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold', 'Color', co(1,:));
end
if mdl3_hi.Coefficients.pValue(2)>pValThresh
    text(.55,.05,sprintf('$p=%.2f, R^2=%.2f$', mdl3_hi.Coefficients.pValue(2), mdl3_hi.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold', 'Color', co(2,:));
else
    text(.55,.05,sprintf('$p<%.4f, R^2=%.2f$', pValThresh, mdl3_hi.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize_txt, 'interpreter', 'latex', 'FontWeight', 'bold', 'Color', co(2,:));
end


ll(1)= plot(nan, nan, 'o-', 'color', co(1,:), 'MarkerSize', mrkSize, 'LineWidth', lw3);
ll(2)= plot(nan, nan, 'd-', 'color', co(2,:), 'MarkerSize', mrkSize, 'LineWidth', lw3);
legend(ll, 'NH', 'HI', 'Location', 'northwest');

mdl_lf_hf= fitlm(lf_x_vals, hf_x_vals);

set(gcf, 'Units', 'inches', 'Position', [1 1 7 9]);

saveas(gcf, 'Figure_Out/all_corr_sfr_limited_2panel', 'png');

saveas(gcf, [LatexDir '/all_corr_sfr_limited_2panel'], 'epsc');