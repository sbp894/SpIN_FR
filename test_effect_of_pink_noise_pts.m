clear;
clc;

pts_file_S= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_10_31-Q366_PTS_SFR/a0002_FFR_SNRenvSSN_Stim_S_P_nType0_atten10.mat';
pts_pinkmask_data_dir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_11_08-Q366_SFR_PTS_pilot2_pink/';
chinID= 366;
stimDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/stim/';
fig_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/maskerSFR/';

if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end


allfiles= dir([pts_pinkmask_data_dir 'a*SFR*.mat']);
nTypes= {'HP', 'LP'};

condStruct.fNames= {allfiles.name}';
condStruct.nType= nTypes(1+cellfun(@(x) contains(x, '_lp_'), {allfiles.name}'))';
snr_start= cellfun(@(x) strfind(x, '_SNR'), {allfiles.name}')+4;
snr_end= cellfun(@(x) strfind(x, '_num'), {allfiles.name}')-1;
condStruct.snr= cell2mat(cellfun(@(x,y,z) str2double(strrep(x(y:z), 'm', '-')), {allfiles.name}', num2cell(snr_start), num2cell(snr_end), 'uniformoutput', false));

uniq_snrs= unique([condStruct.snr]);
% uniq_nType= unique(condStruct.nType); 
uniq_nType= {'HP'}; % ignore LP
saveFigs=1;
% tStart= 0; tEnd= 1.3;
tStart= .3; tEnd= .7;
% tStart= .18; tEnd= .72;
tStart_whole= 0; tEnd_whole= 1.3;

fig_save_dir_subdir= [fig_save_dir sprintf('t%.0fto%.0f_ms/', tStart*1e3, tEnd*1e3)];

raw_power_env= nan(length(uniq_snrs)+1, length(uniq_nType)); % +1 for clean speech
frac_power_env= nan(length(uniq_snrs)+1, length(uniq_nType)); % +1 for clean speech
raw_power_tfs= nan(length(uniq_snrs)+1, length(uniq_nType)); % +1 for clean speech
frac_power_tfs= nan(length(uniq_snrs)+1, length(uniq_nType)); % +1 for clean speech

for snrVar= 1:length(uniq_snrs)
    curSNR= uniq_snrs(snrVar);
    for nTypeVar= 1:length(uniq_nType)
        cur_nType= uniq_nType{nTypeVar};
        
        cur_inds= find(([condStruct.snr]==curSNR) & (strcmp(condStruct.nType, cur_nType)));
        if ~isempty(cur_inds)
            stim_fName= [stimDir strrep(sprintf('Stim_SN_%s_Pink_SNR%d_num1.wav', lower(cur_nType), curSNR), '-', 'm')];
            [sig, fs_sig]= audioread(stim_fName);
            
            sn_data_cell= cell(length(cur_inds), 2);
            nPairs_actual= nan(length(cur_inds), 1);
            for fileVar= 1:length(cur_inds)
                fileInd= cur_inds(fileVar);
                temp_data= load([pts_pinkmask_data_dir allfiles(fileInd).name]);
                temp_data = temp_data.data;
                sn_data_cell{fileVar, 1}= temp_data.AD_Data.AD_Avg_PO_V{1};
                sn_data_cell{fileVar, 2}= temp_data.AD_Data.AD_Avg_NP_V{1};
                
                nPairs_actual(fileVar)= temp_data.Stimuli.RunLevels_params.nPairs_actual;
            end
            
            s_atten=temp_data.Stimuli.atten_dB;
            
            fs_data= temp_data.Stimuli.RPsamprate_Hz;
            t_data= (1:length(sn_data_cell{1,1}))/fs_data;
            inds2load= t_data>tStart_whole & t_data<tEnd_whole;
            sn_data_pos= zeros(1, sum(inds2load));
            sn_data_neg= zeros(1, sum(inds2load));
            
            
            for fileVar=1:length(cur_inds)
                temp_pos= sn_data_cell{fileVar, 1};
                sn_data_pos= sn_data_pos + temp_pos(inds2load)*nPairs_actual(fileVar)/sum(nPairs_actual);
                
                temp_neg= sn_data_cell{fileVar, 2};
                sn_data_neg= sn_data_neg + temp_neg(inds2load)*nPairs_actual(fileVar)/sum(nPairs_actual);
            end
            
            temp_s= load(pts_file_S);
            temp_s= temp_s.data;
            s_data_pos= temp_s.AD_Data.AD_Avg_PO_V{1};
            s_data_neg= temp_s.AD_Data.AD_Avg_NP_V{1};
            
            
            curFilt= get_filter(fs_data);
            
            initialRampDur= 50e-3;
            ramp_nSamples= round(initialRampDur*fs_data);
            sn_rampHamming= hamming(2*ramp_nSamples)';
            sn_rampVector= [sn_rampHamming(1:ramp_nSamples), ones(1, length(sn_data_pos)-length(sn_rampHamming)) sn_rampHamming(ramp_nSamples+1:end)];
            sn_data_pos_filt= filtfilt(curFilt, sn_data_pos.*sn_rampVector).*sn_rampVector;
            sn_data_neg_filt= filtfilt(curFilt, sn_data_neg.*sn_rampVector).*sn_rampVector;
            sn_data_env= (sn_data_pos_filt+sn_data_neg_filt)/2;
            sn_data_tfs= (sn_data_pos_filt-sn_data_neg_filt)/2;
            
            s_rampHamming= hamming(2*ramp_nSamples)';
            s_rampVector= [s_rampHamming(1:ramp_nSamples), ones(1, length(s_data_pos)-length(s_rampHamming)) s_rampHamming(ramp_nSamples+1:end)];
            s_data_pos_filt= filtfilt(curFilt, s_data_pos.*s_rampVector).*s_rampVector;
            s_data_neg_filt= filtfilt(curFilt, s_data_neg.*s_rampVector).*s_rampVector;
            s_data_env= (s_data_pos_filt+s_data_neg_filt)/2;
            s_data_tfs= (s_data_pos_filt-s_data_neg_filt)/2;
            
            figHan= 1;
            fName= strrep(sprintf('Q%d_pts_SNR%d_%s_s_vs_sn', chinID, curSNR, cur_nType), '-', 'm');
            ttlStr= sprintf('Q%d, PTS, SNR %d, %s', chinID, curSNR, cur_nType);
            %             create_panel_plot_include_hilbert(figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, sn_data_env, sn_data_tfs, tStart, tEnd, ...
            %                 [fig_save_dir sprintf('t%.0fto%.0f_ms/', tStart*1e3, tEnd*1e3)], fName, ttlStr, saveFigs);
            PSD_struct= create_panel_plot_s_vs_sn(figHan, fs_sig, sig, fs_data, sn_data_pos_filt, sn_data_neg_filt, s_data_pos_filt, s_data_neg_filt, tStart, tEnd, ...
                fig_save_dir_subdir, fName, ttlStr, saveFigs);
            
            raw_power_env(snrVar, nTypeVar)= PSD_struct.raw.ENV.SN;
            raw_power_tfs(snrVar, nTypeVar)= PSD_struct.raw.TFS.SN;
            frac_power_env(snrVar, nTypeVar)= PSD_struct.frac.ENV.SN;
            frac_power_tfs(snrVar, nTypeVar)= PSD_struct.frac.TFS.SN;
            
            if snrVar==length(uniq_snrs)
                raw_power_env(length(uniq_snrs)+1, nTypeVar)= PSD_struct.raw.ENV.S;
                raw_power_tfs(length(uniq_snrs)+1, nTypeVar)= PSD_struct.raw.TFS.S;
                frac_power_env(length(uniq_snrs)+1, nTypeVar)= PSD_struct.frac.ENV.S;
                frac_power_tfs(length(uniq_snrs)+1, nTypeVar)= PSD_struct.frac.TFS.S;          
            end
            
        end
    end
end

%%
raw_power_env= raw_power_env(:,1);
raw_power_tfs= raw_power_tfs(:,1);
frac_power_env= frac_power_env(:,1);
frac_power_tfs= frac_power_tfs(:,1);

figure(11);
clf;
mrkSize= 20;
lw= 4;
fSize= 20;
xlabel_str= cellfun(@(x) num2str(x), num2cell([uniq_snrs; inf]), 'uniformoutput', false);

subplot(211);
hold on;
plot(1:size(raw_power_env,1), raw_power_env, 'v', 'markersize', mrkSize, 'linew', lw);
plot(1:size(raw_power_tfs,1), raw_power_tfs, '^', 'markersize', mrkSize, 'linew', lw);
% legend('ENV', 'TFS');
grid on;
ylabel('$RAW_{power} (dB)$', 'interpreter', 'latex');
set(gca,'fontsize', fSize, 'xtick', 1:size(raw_power_env,1), 'xticklabel', xlabel_str);
title('PTS- Q366');

subplot(212);
hold on;
plot(1:size(frac_power_env,1), frac_power_env, 'v', 'markersize', mrkSize, 'linew', lw);
plot(1:size(frac_power_tfs,1), frac_power_tfs, '^', 'markersize', mrkSize, 'linew', lw);
legend('ENV', 'TFS', 'location', 'southeast');
grid on;
ylabel('$FRAC_{power} (dB)$', 'interpreter', 'latex');
set(gca,'fontsize', fSize, 'xtick', 1:size(raw_power_env,1), 'xticklabel', xlabel_str);
xlabel('HP Masking SNR (dB)');
set(gcf, 'units', 'inches', 'position', [1 1 10 6]);

fName_summary= sprintf('Q%d_pts_pink_summary_%s', chinID, 'HP');
saveas(gcf, [fig_save_dir_subdir fName_summary], 'png');

function curFilt= get_filter(fs_data)
N_bp_half= 4;
HalfPowerFrequency1=0.5;
HalfPowerFrequency2=1e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end