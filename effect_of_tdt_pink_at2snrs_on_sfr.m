clear;
clc;
clf;

chinID= 371;
data_dir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Data_Out/Artifact_Removed_FFR/SP-2019_05_08-Q371_SFRpink500Hz_NH/';
remove_artifact_here= 0; % because using cleaned data => already artifact removed
if ~remove_artifact_here
    fprintf('Assuming using artifact removed Data\n');
end
tStart= 0; tEnd= 1.3; 

allfiles= dir([data_dir 'a*SFR*.mat']);
allfiles= allfiles(~(contains({allfiles.name}, 'latency') | contains({allfiles.name}, 'artifact')));

all_snrs= cell2mat(cellfun(@(x) str2double(strrep(x(regexp(x, 'snr_')+4 : regexp(x, '_atn')-1), 'm', '-')), {allfiles.name}, 'uniformoutput', false));
all_snrs(isnan(all_snrs))= [];
all_snrs= fliplr(unique(all_snrs));

stim_fName= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2017_11_02-Q325_AN_NH/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav';
[sig, fs_sig]= audioread(stim_fName);
snrVar= 1;
[s_data_pos_filt, s_data_neg_filt, s_nf_pos_filt, s_nf_neg_filt, fs_data]= tdt_pink_helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);

snrVar= find(all_snrs==0);
[sn_data_pos_filt_p10, sn_data_neg_filt_p10, sn_nf_pos_filt_p10, sn_nf_neg_filt_p10, fs_data_p10]= tdt_pink_helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);

snrVar= find(all_snrs==-20);
[sn_data_pos_filt_m20, sn_data_neg_filt_m20, sn_nf_pos_filt_m20, sn_nf_neg_filt_m20, fs_data_m20]= tdt_pink_helper.get_filtered_ffr(snrVar, allfiles, data_dir, all_snrs, remove_artifact_here);

sig= gen_resample(sig, fs_sig, fs_data_p10);
fs_sig= fs_data_p10;

plot_S_sig= 1;
legHan1= create_panel_plot_s_vs_sn_tdt_at2snrs(1, fs_sig, sig, fs_data_p10, sn_data_pos_filt_p10, sn_data_neg_filt_p10, s_data_pos_filt, s_data_neg_filt, ...
    tStart, tEnd, plot_S_sig);

plot_S_sig= 0;
legHan2= create_panel_plot_s_vs_sn_tdt_at2snrs(1, fs_sig, sig, fs_data_m20, sn_data_pos_filt_m20, sn_data_neg_filt_m20, s_data_pos_filt, s_data_neg_filt, ...
    tStart, tEnd, plot_S_sig);

legHans= [legHan1; legHan2];
lg= legend(legHans, 'Stim', 'Quiet', 'SNR 0 dB', 'SNR -20 dB', 'box', 'off', 'location', 'southwest');
lg.FontSize= 14;
set(gcf, 'units', 'centimeters', 'position', [50 3 30 15], 'Renderer','painters');

LatexDir= '/home/parida/Dropbox/Articles/Loss_of_tonotopy_in_HI_FFR/figures/';
fName_summary= sprintf('nh_tdt_pink_example');
saveas(gcf, [LatexDir fName_summary], 'png');
saveas(gcf, [LatexDir fName_summary], 'tiff');