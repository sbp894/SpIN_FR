clear;
clc;

AllChinIDs= [358 360 366 367 369 370 371 373 374 379];
% AllChinIDs= [369 373];
% AllChinIDs= 369; %[373 374];
hearingTypes= {'NH', 'PTS'};
% AllChinIDs= 369;


CodesDirs= {'/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/helper', ...
    '/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/continuous'};
rmpath(CodesDirs{:});
% findpeaks code in the above folder shadows matlab's findpeaks


latexDir= '/home/parida/Dropbox/Articles/Loss_of_tonotopy_in_HI_FFR/figures/';
saveAgain=1;
saveLatex= 0;

fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/SFR_corr_latency/');
if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end

data_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/SFR_corr_latency/');
if ~isfolder(data_save_dir)
    mkdir(data_save_dir);
end

fixed_acoustic_delay= 2.39e-3;
filt_ffr.BPlow=100;
filt_ffr.BPhigh=.5e3;
filt_sig.LPco=filt_ffr.BPhigh;
filt_sig.BPlow=.5e3;
filt_sig.BPhigh=8e3;

[sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');
t_sig= (1:length(sig))/fs;
stim_dur= length(sig)/fs;


windowLength= 64e-3;
fracOverLap= .75;
NW= 1.5; % NW= dur * f_res => f_res= NW/dur. If NW=1.5, dur=50ms, f_res= 30 Hz

fracSlide= 1-fracOverLap;
tSlide= fracSlide*windowLength;
nSegs= 1 + floor((stim_dur-windowLength)/(tSlide));

dpss_yRange= 100;
plot_dpss_each_iter= false;
lw= 2;
lw2= 3;
nSProws=2;
fSize= 20;
nSPcols=1;

t_latency= 5e-3;

lat_env_nh= [];
lat_tfs_nh= [];
lat_env_hi= [];
lat_tfs_hi= [];

for chinVar= 1:length(AllChinIDs)
    restart_fig= 1;
    
    
    chinID= AllChinIDs(chinVar);
    
    if chinID==191
        warning('Don''t have data in response to speech, Collected data to N twice. Recollect data');
        pause
    end
    
    RootDataDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Artifact_Removed_FFR/';
    
    allFiles= dir([RootDataDir '*' num2str(chinID) '*SFR*']);
    
    for typeVar= 1:length(hearingTypes)
        postFix= hearingTypes{typeVar};
        validInds= contains({allFiles.name}, postFix) & ~contains({allFiles.name}, 'pink');
        
        figure(1);
        clf;
        co= get(gca, 'colororder');
        if chinID==369 && strcmp(postFix,'PTS')
            validInds= 0;
        end
        
        if any(validInds)
            data_dir= [RootDataDir allFiles(validInds).name filesep];
            fprintf('Working on %s\n', data_dir);
            s_files= dir([data_dir 'a*_S_*']);
            %% clean speech
            s_data_cell= cell(length(s_files), 2);
            nPairs_actual= nan(length(s_files), 1);
            for sfile_var=1:length(s_files)
                temp_data= load([data_dir s_files(sfile_var).name]);
                temp_data = temp_data.data;
                s_data_cell{sfile_var, 1}= temp_data.AD_Data.AD_Avg_NP_V{1}; % Note PO and NP were switched in all FFR before 22 July 2019
                s_data_cell{sfile_var, 2}= temp_data.AD_Data.AD_Avg_PO_V{1};
                
                nPairs_actual(sfile_var)= temp_data.Stimuli.RunLevels_params.nPairs_actual;
            end
            
            s_atten=temp_data.Stimuli.atten_dB;
            
            [s_sig, fs_sig]= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_09_11-Q362_AN_PTS/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');
            t_sig= (1:length(s_sig))/fs_sig;
            
            s_data_pos= zeros(1, length(s_data_cell{sfile_var,1}));
            s_data_neg= zeros(1, length(s_data_cell{sfile_var,2}));
            fs_data= temp_data.Stimuli.RPsamprate_Hz;
            
            
            for i=1:length(s_files)
                s_data_pos= s_data_pos + s_data_cell{i, 1}*nPairs_actual(i)/sum(nPairs_actual);
                s_data_neg= s_data_neg + s_data_cell{i, 2}*nPairs_actual(i)/sum(nPairs_actual);
            end
            
            initialRampDur= 20e-3;
            ramp_nSamples= round(initialRampDur*fs_data);
            rampHamming= hamming(2*ramp_nSamples)';
            rampVector= [rampHamming(1:ramp_nSamples), ones(1, length(s_data_pos)-length(rampHamming)) rampHamming(ramp_nSamples+1:end)];
            bpFilt_ffr= get_bp_filter(fs_data, filt_ffr.BPlow, filt_ffr.BPhigh);
            
            s_data_pos_filt= filtfilt(bpFilt_ffr, s_data_pos.*rampVector).*rampVector;
            s_data_neg_filt= filtfilt(bpFilt_ffr, s_data_neg.*rampVector).*rampVector;
            s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
            
            if ismember(chinID, [373 374 379])
                s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;
            else
                s_data_tfs= (s_data_neg_filt-s_data_pos_filt)/2;
            end
            t_data= (1:length(s_data_env))/fs_data;
            
            sig_rs= gen_resample(sig, fs_sig, fs_data);
            t_sig_rs= (1:length(sig_rs))/fs_data;
            
            bpFilt_sig= get_bp_filter(fs_data, filt_sig.BPlow, filt_sig.BPhigh);
            lpFilt_sig= get_lp_filter(fs_data, filt_sig.LPco);
            
            sig_rs_env= filtfilt(lpFilt_sig, abs(hilbert(filtfilt(bpFilt_sig, sig_rs))));
            sig_rs_tfs= filtfilt(lpFilt_sig, sig_rs);
            sig_rs_tfs(sig_rs_tfs<0)= 0; % HWR
            
            s_data_env_trim= s_data_env(1:round(stim_dur*fs_data));
            s_data_tfs_trim= s_data_tfs(1:round(stim_dur*fs_data));
            
            [ccf_Renv_Senv, ~]= xcorr(s_data_env_trim, sig_rs_env); % stimulus should be the second argIN
            [ccf_Rtfs_Stfs, delay]= xcorr(s_data_tfs_trim, sig_rs_tfs); % stimulus should be the second argIN
%             ccf_Rtfs_Stfs= abs(ccf_Rtfs_Stfs);
%             warning('debugging');
            
            delay= delay/fs_data;
            delay= delay-fixed_acoustic_delay;
            
            minDelay= 2e-3;
            MinPeakWidth= round(1e-3*fs_data);
            [env_peak_val, env_peak_ind]= findpeaks(ccf_Renv_Senv, 'MinPeakWidth', MinPeakWidth);
            env_cand_ind= find((delay(env_peak_ind)>minDelay), 1);
            env_peak_x= delay(env_peak_ind(env_cand_ind))*1e3;
            env_peak_y= env_peak_val(env_cand_ind);
            
            [tfs_peak_val, tfs_peak_ind]= findpeaks(ccf_Rtfs_Stfs, 'MinPeakWidth', MinPeakWidth);
            tfs_cand_ind= find((delay(tfs_peak_ind)>minDelay), 1);
            tfs_peak_x= delay(tfs_peak_ind(tfs_cand_ind))*1e3;
            tfs_peak_y= tfs_peak_val(tfs_cand_ind);
            
            
            figure(1);
            clf;
            hold on;
            lw= 2;
            lw2= 3;
            fSize= 20;
            mrkSize= 16;
            plot(delay*1e3, ccf_Renv_Senv, '-', delay*1e3, ccf_Rtfs_Stfs, '-', 'linew', lw);
            plot(env_peak_x, env_peak_y, 'v', 'color', co(1,:), 'markersize', mrkSize, 'linew', lw2);
            plot(tfs_peak_x, tfs_peak_y, 'v', 'color', co(2,:), 'markersize', mrkSize, 'linew', lw2);
            text(env_peak_x, env_peak_y*1.15, sprintf('$delay_{ENV}=%.1f ms$', env_peak_x), 'fontsize', fSize, 'interpreter', 'latex', 'color', co(1,:));
            text(tfs_peak_x, tfs_peak_y*1.15, sprintf('$delay_{TFS}=%.1f ms$', tfs_peak_x), 'fontsize', fSize, 'interpreter', 'latex', 'color', co(2,:));
            
            if strcmp(postFix, 'NH')
                lat_env_nh= [lat_env_nh; [round(chinID), env_peak_x]]; %#ok<*AGROW>
                lat_tfs_nh= [lat_tfs_nh; [round(chinID), tfs_peak_x]];
                
            elseif strcmp(postFix, 'PTS')
                lat_env_hi=[lat_env_hi; [round(chinID), env_peak_x]];
                lat_tfs_hi= [lat_tfs_hi; [round(chinID), tfs_peak_x]];
            end
            
            grid off;
            
            legend('ENV', 'TFS', 'box', 'off');
            
            title(sprintf('Q%d | %s', chinID, postFix));
            xlabel('Delay (ms)');
%             ylabel('CCF');
            ylabel({'Cross-correlation'; 'Function'});
            xlim([-5 25]);
            set(gca, 'fontsize', fSize);
            set(gcf, 'units', 'inches', 'position', [1 1 10 5]);
            
            ylim_new= 1.3*max(abs([ccf_Renv_Senv(:); ccf_Rtfs_Stfs(:)]));
            ylim([-.85*ylim_new ylim_new]);
            
            fName= sprintf('%sQ%d_%s_latency', fig_save_dir, chinID, postFix);
            if saveAgain
                saveas(gcf, fName, 'png');
            end
            fName= sprintf('%sQ%d_%s_latency', latexDir, chinID, postFix);
            
            if saveLatex & contains(fName, 'Q369_NH_latency')
                title('')
                saveas(gcf, fName, 'epsc');
            end
            figure(11);
            clf
            plot(t_sig_rs, sig_rs_tfs)
            yyaxis right
            plot(t_data-fixed_acoustic_delay, s_data_tfs)
            
            title(num2str(chinID))
            
        end
    end
end

save([data_save_dir 'latency_data.mat'], 'lat_env_nh', 'lat_tfs_nh', 'lat_env_hi', 'lat_tfs_hi');

function bpFilt= get_bp_filter(fs_data, HalfPowerFrequency1, HalfPowerFrequency2)
N_bp_half= 4;


bpFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end

function lpFilt= get_lp_filter(fs_data, PassbandFrequency)
N_lp_half= 4;
PassbandRipple= .2;

lpFilt = designfilt('lowpassiir','FilterOrder', N_lp_half, ...
    'PassbandFrequency', PassbandFrequency, 'PassbandRipple', PassbandRipple, ...
    'SampleRate',fs_data);

end