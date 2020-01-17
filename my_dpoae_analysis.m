function out_DPOAE_data= my_dpoae_analysis(fName, plotVar)

if ~exist('plotVar', 'var')
    plotVar= false;
end

run(fName);
dpData= ans;

% clear curData
out_DPOAE_data= repmat(struct('freq1', nan, 'freq2', nan, 'freq_dp', nan, 'dp_amp', nan), size(dpData.DpoaeData, 1), 1);

newValue= num2cell(dpData.DpoaeData(:,1));
[out_DPOAE_data.freq1]= newValue{:};

newValue= num2cell(dpData.DpoaeData(:,2));
[out_DPOAE_data.freq2]= newValue{:};

newValue= num2cell(dpData.DpoaeData(:,3));
[out_DPOAE_data.freq_dp]= newValue{:};

psd_freq= dpData.Dpoaefreqs;
neighborhood_freq_inds= 15;

for freqVar= 1:size(dpData.DpoaeData, 1)
    curPSD= dpData.DpoaeSpectra(freqVar, :);
    cur_dpFreq= out_DPOAE_data(freqVar).freq_dp;
    
    
    
    ind_dp= dsearchn(psd_freq', cur_dpFreq);
    inds_lin_fit= ind_dp-neighborhood_freq_inds:ind_dp+neighborhood_freq_inds;
    mad_val= mad(curPSD(inds_lin_fit));
    
    mdl= fitlm(psd_freq(inds_lin_fit), curPSD(inds_lin_fit));
    fitted_data= predict(mdl,psd_freq(inds_lin_fit)');
    
    [~, max_ind_rel]= max(curPSD(inds_lin_fit));
    max_ind= ind_dp-neighborhood_freq_inds-1+max_ind_rel;
    
    
    if plotVar
        figure(111);
        clf;
        
        plot(psd_freq, curPSD);
        set(gca, 'xscale', 'log');
        ylim_cur= ylim;
        hold on;
        plot(psd_freq(inds_lin_fit), fitted_data+mad_val);
        plot(repmat(cur_dpFreq, 1, 2), ylim_cur, 'k');
        xlim([cur_dpFreq*.9 out_DPOAE_data(freqVar).freq2*1.1]);
    end
    
    if curPSD(max_ind)> (predict(mdl,psd_freq(max_ind))+mad_val)
        out_DPOAE_data(freqVar).dp_amp= curPSD(max_ind);
        if plotVar
            plot(psd_freq(max_ind), curPSD(max_ind), 'r*', 'markersize', 18);
        end
    end
end