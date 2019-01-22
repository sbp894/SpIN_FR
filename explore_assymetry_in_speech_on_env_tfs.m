clear;
clc;

all_f_highs= [.5e3 1e3 2e3 4e3 8e3];

for freqVar=1:length(all_f_highs)
    f_high= all_f_highs(freqVar);
    [sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');
    t_sig= (1:length(sig))/fs;
    
    fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/explore_assymetry/');
    if ~isdir(fig_save_dir)
        mkdir(fig_save_dir);
    end
    
    bw_ratio= logspace(log10(.01), log10(.99), 12);
    plot_stages= 1;
    
    N_bp_half=4;
    

    
    env_2_tfs_ratio= nan(length(bw_ratio), 1);
    bw_in_Hz= nan(length(bw_ratio), 1);
    
    for bwVar= 1:length(bw_ratio)
        cur_bw_ratio= bw_ratio(bwVar);
        f_low= (1-cur_bw_ratio)*f_high;
    
        HalfPowerFrequency2= min(.5, .9*f_low);
        lpFilt = designfilt('bandpassiir','FilterOrder',40, ...get(p,props)
        'HalfPowerFrequency1',1,'HalfPowerFrequency2',.05e3, ...
        'SampleRate',fs);
    
        cur_bp_Filt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
            'HalfPowerFrequency1',f_low,'HalfPowerFrequency2',f_high, ...
            'SampleRate',fs);
        bp_out_pos= filter(cur_bp_Filt, sig);
        hwr_out_pos= bp_out_pos .* double(bp_out_pos>0);
        env_out_pos= filter(lpFilt, hwr_out_pos);
        
        bp_out_neg= filter(cur_bp_Filt, -sig);
        hwr_out_neg= bp_out_neg .* double(bp_out_neg>0);
        env_out_neg= filter(lpFilt, hwr_out_neg);
        
        if plot_stages
            nSProws=4;
            nSPcols=1;
            
            figure(1);
            ax(1)= subplot(nSProws, nSPcols, 1);
            plot(t_sig, sig);
            
            ax(2)= subplot(nSProws, nSPcols, 2);
            plot(t_sig, bp_out_pos, t_sig, bp_out_neg);
            
            ax(3)= subplot(nSProws, nSPcols, 3);
            plot(t_sig, hwr_out_pos, t_sig, hwr_out_neg);
            
            ax(4)= subplot(nSProws, nSPcols, 4);
            plot(t_sig, env_out_pos, t_sig, env_out_neg);
            
            linkaxes(ax, 'x');
        end
        env_out= (env_out_pos+env_out_neg)/2;
        tfs_out= (env_out_pos-env_out_neg)/2;
        env_2_tfs_ratio(bwVar)= rms(env_out)/rms(tfs_out);
        bw_in_Hz(bwVar)= f_high-f_low;
    end
    
    %%
    figure(2);
    lw=4;
    mrkSize=16;
    fSize= 20;
    loglog(bw_in_Hz, env_2_tfs_ratio, '-d', 'linew', lw, 'markersize', mrkSize);
    xlabel('BW (Hz)');
    grid on;
    ylabel('audio rms(ENV)/rms(TFS) ');
    title(sprintf('measure of assymetry (HWR+LPF), f-up-co=%.1fk', f_high/1e3));
    set(gca, 'fontsize', fSize);
    
    set(2, 'units', 'normalized', 'position', [.1 .1 .8 .8]);
    figName= sprintf('assymetry_upper_freq_co_%.0f', f_high);
    saveas(gcf, [fig_save_dir figName], 'tiff');
end
fprintf('Conclusion: if lowest frequency in filter is above HWR cut-off, \n envelopes identical for both polarities. if not, TFS passes through \n and assymetry in response gets picked up. \n');
