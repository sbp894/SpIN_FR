clear;
close all;
clc;


[sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');
outputDir.stim= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/stim/';
if ~isdir(outputDir.stim)
    mkdir(outputDir.stim);
end
outputDir.psd= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/psd/';
if ~isdir(outputDir.psd)
    mkdir(outputDir.psd);
end

sig_len= length(sig);
sig_dur= sig_len/fs;
t_data= (1:sig_len)/fs;
osSpl= 20*log10(rms(sig)/(20e-6));

cn= dsp.ColoredNoise('Color','pink', 'SamplesPerFrame', sig_len);
pn1= cn();


%%
% figure(1); clf;
% subplot(211);
% plot(t_data, pn1(1:sig_len));
% 
% subplot(212);
% plot_dpss_psd(pn1(1:sig_len), fs, 'NW', 10);

%%
% going to do 16 reps of 16 different noise iterations
num_of_NoisySig= 16;
all_snr= -20:10:20;

N_HP= 40;
F3dB_HP= 1.5e3;
Ast_HP= 100;
hd_HP= fdesign.highpass('N,F3dB,Ast', N_HP, F3dB_HP, Ast_HP, fs);
fmethds_HP= designmethods(hd_HP);
hpFilt= design(hd_HP, fmethds_HP{1});

N_LP= 40;
F3dB_LP= 1.5e3;
Ast_LP= 100;
hd_LP= fdesign.lowpass('N,F3dB,Ast', N_LP, F3dB_LP, Ast_LP, fs);
fmethds_LP= designmethods(hd_LP);
lpFilt= design(hd_LP, fmethds_LP{1});


fSize= 16;
for snrVar= 1:length(all_snr)
    cur_snr= all_snr(snrVar);
    
    for nIter= 1: num_of_NoisySig
        [~, noise_pink]=create_noisy_signal(sig, cur_snr, 'pink');
        
        hp_noise= filter(hpFilt, noise_pink);
        sn_pink_hp= create_noisy_signal(sig, cur_snr, hp_noise);
        
        %         sound(sn_pink_hp, fs);
        %         pause(length(sn_pink_hp)/fs);
        
        lp_noise= filter(lpFilt, noise_pink);
        sn_pink_lp= create_noisy_signal(sig, cur_snr, lp_noise);
        
        do_plotting(1, sig, fs, noise_pink, hp_noise, sn_pink_hp, lp_noise, sn_pink_lp, fSize);
        %         sound(sn_pink_lp, fs);
        %         pause(length(sn_pink_lp)/fs);
        
        fName_hp= strrep(sprintf('Stim_SN_hp_Pink_SNR%d_num%d.wav', cur_snr, nIter), '-', 'm');
        dummy_aud_wrt([outputDir.stim fName_hp], sn_pink_hp, fs);
        
        fName_lp= strrep(sprintf('Stim_SN_lp_Pink_SNR%d_num%d.wav', cur_snr, nIter), '-', 'm');
        audiowrite([outputDir.stim fName_lp], sn_pink_lp, fs);
        
        set(gcf, 'units', 'normalized', 'position', [.1 .1 .8 .8]);
        figName= strrep(sprintf('Stim_SN_Pink_SNR%d_num%d', cur_snr, nIter), '-', 'm');
        dummy_saveas(1, [outputDir.psd figName]);
    end
end

function dummy_aud_wrt(fName, sig, fs)
audiowrite(fName, sig, fs);
end

function dummy_saveas(figHan, figName)
saveas(figHan, figName);
saveas(figHan, figName, 'tiff');
end

function do_plotting(figHan, sig, fs, noise_pink, hp_noise, sn_pink_hp, lp_noise, sn_pink_lp, fSize)
figure(figHan);
clf;
subplot(121);
hold on;
plot_dpss_psd(sig, fs, 'NW', 4);
plot_dpss_psd(sn_pink_hp, fs, 'NW', 4);
plot_dpss_psd(noise_pink, fs, 'NW', 4);
plot_dpss_psd(hp_noise, fs, 'NW', 4);
legend('S', 'SN', 'N', 'N-HP');
set(gca, 'fontsize', fSize);

subplot(122);
hold on;
plot_dpss_psd(sig, fs, 'NW', 4);
plot_dpss_psd(sn_pink_lp, fs, 'NW', 4);
plot_dpss_psd(noise_pink, fs, 'NW', 4);
plot_dpss_psd(lp_noise, fs, 'NW', 4);
legend('S', 'SN', 'N', 'N-LP');
set(gca, 'fontsize', fSize);
end