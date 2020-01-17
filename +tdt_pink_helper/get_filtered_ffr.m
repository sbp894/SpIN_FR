function [sn_data_pos_filt, sn_data_neg_filt, nf_data_pos_filt, nf_data_neg_filt, fs_new]= get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact)
curSNR= all_snrs(snrVar);
curFileInd= contains({allfiles.name}, strrep(['snr_' num2str(curSNR)], '-', 'm'));
curFile= allfiles(curFileInd).name;

fprintf('Using %s for %d dB SNR\n', curFile, curSNR);


temp= load([data_dir curFile]);

fs_org= temp.data.Stimuli.RPsamprate_Hz;
fs_new= 4e3;

if remove_artifact
    fMax= 3e3;
    plotVar= 0;
    sn_data_pos= remove_artifact_ffr(temp.data.AD_Data.AD_Avg_PO_V{1}, fs_org, plotVar, fMax);
    sn_data_neg= remove_artifact_ffr(temp.data.AD_Data.AD_Avg_NP_V{1}, fs_org, plotVar, fMax);
    nf_data_pos= remove_artifact_ffr(temp.data.AD_Data.AD_NF_PO_V{1}, fs_org, plotVar, fMax);
    nf_data_neg= remove_artifact_ffr(temp.data.AD_Data.AD_NF_NP_V{1}, fs_org, plotVar, fMax);
else
    sn_data_pos= temp.data.AD_Data.AD_Avg_PO_V{1};
    sn_data_neg= temp.data.AD_Data.AD_Avg_NP_V{1};
    nf_data_pos= temp.data.AD_Data.AD_NF_PO_V{1};
    nf_data_neg= temp.data.AD_Data.AD_NF_NP_V{1};
end


initialRampDur= 50e-3;
ramp_nSamples= round(initialRampDur*fs_org);
bpFilt= tdt_pink_helper.get_filter(fs_org);

rampHamming= hamming(2*ramp_nSamples)';
rampVector= [rampHamming(1:ramp_nSamples), ones(1, length(sn_data_pos)-length(rampHamming)) rampHamming(ramp_nSamples+1:end)];

%
sn_data_pos_filt= filtfilt(bpFilt, sn_data_pos.*rampVector);
sn_data_neg_filt= filtfilt(bpFilt, sn_data_neg.*rampVector);
nf_data_pos_filt= filtfilt(bpFilt, nf_data_pos.*rampVector);
nf_data_neg_filt= filtfilt(bpFilt, nf_data_neg.*rampVector);

%
sn_data_pos_filt= gen_resample(sn_data_pos_filt, fs_org, fs_new);
sn_data_neg_filt= gen_resample(sn_data_neg_filt, fs_org, fs_new);
nf_data_pos_filt= gen_resample(nf_data_pos_filt, fs_org, fs_new);
nf_data_neg_filt= gen_resample(nf_data_neg_filt, fs_org, fs_new);

end