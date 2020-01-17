clear;

fName_nh= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2019_03_19-Q371_SFR_pink_NH/a0008_SFR_pink_S_snr_120_atn10_latency.mat';
latency_nh= check_acosutic_delay(fName_nh)

fName_hi= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2019_03_19-Q370_SFR_pink_PTS/a0008_SFR_pink_S_snr_120_atn10_latency.mat';
latency_hi= check_acosutic_delay(fName_hi)

avg_latency = (latency_nh.pos + latency_nh.neg + latency_hi.pos + latency_hi.neg)/4
