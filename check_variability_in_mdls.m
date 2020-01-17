function check_variability_in_mdls(pool_ratio_lf_to_hf_audio, pool_ratio_tfs_to_env_ffr, pool_lf_power_audio, pool_hf_power_audio, pool_env_power_ffr, pool_tfs_power_ffr, nhInds, hiInds)


nh_slopes= struct('ratio', nan(1, length(nhInds)), 'tfs', nan(1, length(nhInds)), 'env', nan(1, length(nhInds)));
for chinVar= 1:length(nhInds)
    curInd= nhInds(chinVar);
    
    % ratio 
    x= pool_ratio_lf_to_hf_audio(curInd, :);
    y= pool_ratio_tfs_to_env_ffr(curInd, :);
    mdl1= fitlm(x, y);
    nh_slopes.ratio(chinVar)= mdl1.Coefficients.Estimate(2);
    
    % tfs 
    x= dbspl(sqrt(pool_lf_power_audio(curInd, :)));
    y= db(sqrt(pool_tfs_power_ffr(curInd, :)));
    mdl2= fitlm(x, y);
    nh_slopes.tfs(chinVar)= mdl2.Coefficients.Estimate(2);
    
    % env    
    x= dbspl(sqrt(pool_hf_power_audio(curInd, :)));
    y= db(sqrt(pool_env_power_ffr(curInd, :)));
    mdl3= fitlm(x, y);
    nh_slopes.env(chinVar)= mdl3.Coefficients.Estimate(2);
end

hi_slopes= struct('ratio', nan(1, length(hiInds)), 'tfs', nan(1, length(hiInds)), 'env', nan(1, length(hiInds)));
for chinVar= 1:length(hiInds)
    curInd= hiInds(chinVar);
    
    % ratio 
    x= pool_ratio_lf_to_hf_audio(curInd, :);
    y= pool_ratio_tfs_to_env_ffr(curInd, :);
    mdl1= fitlm(x, y);
    hi_slopes.ratio(chinVar)= mdl1.Coefficients.Estimate(2);
    
    % tfs 
    x= dbspl(sqrt(pool_lf_power_audio(curInd, :)));
    y= db(sqrt(pool_tfs_power_ffr(curInd, :)));
    mdl2= fitlm(x, y);
    hi_slopes.tfs(chinVar)= mdl2.Coefficients.Estimate(2);
    
    % env    
    x= dbspl(sqrt(pool_hf_power_audio(curInd, :)));
    y= db(sqrt(pool_env_power_ffr(curInd, :)));
    mdl3= fitlm(x, y);
    hi_slopes.env(chinVar)= mdl3.Coefficients.Estimate(2);
end

fprintf('---------NH---------\n ratio | tfs | env \n');
[std(nh_slopes.ratio) std(nh_slopes.tfs) std(nh_slopes.env)]
fprintf('---------HI---------\n');
[std(hi_slopes.ratio) std(hi_slopes.tfs) std(hi_slopes.env)]