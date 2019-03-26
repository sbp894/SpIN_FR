clear;
clc;

pts_data_dir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_11_08-Q366_SFR_PTS_pilot2/';
chinID= 366;
stimDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/stim/';
fig_save_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/maskerSFR/';

if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end


allfiles= dir([pts_data_dir 'a*SFR*.mat']);
nTypes= {'HP', 'LP'};

condStruct.fNames= {allfiles.name}';
condStruct.nType= nTypes(1+cellfun(@(x) contains(x, '_lp_'), {allfiles.name}'))';
snr_start= cellfun(@(x) strfind(x, '_SNR'), {allfiles.name}')+4;
snr_end= cellfun(@(x) strfind(x, '_num'), {allfiles.name}')-1;
condStruct.snr= cell2mat(cellfun(@(x,y,z) str2double(strrep(x(y:z), 'm', '-')), {allfiles.name}', num2cell(snr_start), num2cell(snr_end), 'uniformoutput', false));

uniq_snrs= unique([condStruct.snr]);
uniq_nType= unique(condStruct.nType);
saveFigs=1;
% tStart= 0; tEnd= 1.3;
% tStart= .18; tEnd= .72;
tStart= .23; tEnd= .35;
tStart_whole= 0; tEnd_whole= 1.3;

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
                temp_data= load([pts_data_dir allfiles(fileInd).name]);
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
            
            
            curFilt= get_filter(fs_data);
            
            initialRampDur= 50e-3;
            ramp_nSamples= round(initialRampDur*fs_data);
            rampHamming= hamming(2*ramp_nSamples)';
            rampVector= [rampHamming(1:ramp_nSamples), ones(1, length(sn_data_pos)-length(rampHamming)) rampHamming(ramp_nSamples+1:end)];
            data_pos_filt= filtfilt(curFilt, sn_data_pos.*rampVector).*rampVector;
            data_neg_filt= filtfilt(curFilt, sn_data_neg.*rampVector).*rampVector;
            data_env= (data_pos_filt+data_neg_filt)/2;
            data_tfs= (data_pos_filt-data_neg_filt)/2;
            
            
            figHan= 1;
            fName= strrep(sprintf('Q%d_pts_SNR%d_%s', chinID, curSNR, cur_nType), '-', 'm');
            ttlStr= sprintf('Q%d, PTS, SNR %d, %s', chinID, curSNR, cur_nType);
            create_panel_plot_include_hilbert(figHan, fs_sig, sig, fs_data, data_pos_filt, data_neg_filt, data_env, data_tfs, tStart, tEnd, ...
                [fig_save_dir sprintf('t%.0fto%.0f_ms/', tStart*1e3, tEnd*1e3)], fName, ttlStr, saveFigs);
        end
    end
end


function curFilt= get_filter(fs_data)
N_bp_half= 4;
HalfPowerFrequency1=0.5;
HalfPowerFrequency2=1e3;

curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end