clear;
clc;

pRef= 20e-6;
[sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');
nfft= 2^nextpow2(numel(sig));
OAspl= 20*log10(rms(sig)/pRef);
fprintf('Long-term intensity = %.1f dB SPL\n', OAspl);

[Pxx_dB, freq]= plot_dpss_psd(sig, fs, 'nfft', nfft, 'plot', false);

fBoundary= 500; 
inds_low= freq<fBoundary;

Pxx_lin= (10.^(Pxx_dB/10));
var_spectral= sum(Pxx_lin)/nfft*fs;
amp_spectral= sqrt(var_spectral);
OAspl_spectral= 20*log10(amp_spectral/pRef);

Plow_lin= Pxx_lin(inds_low);
Phigh_lin= Pxx_lin(~inds_low);

fprintf('PSD estimated intensity = %.1f dB SPL\n', OAspl_spectral);
fprintf('PSD estimated intensity for LF= %.1f dB SPL\n', 20*log10(sqrt(sum(Plow_lin)/nfft*fs)/pRef));
fprintf('PSD estimated intensity for HF= %.1f dB SPL\n', 20*log10(sqrt(sum(Phigh_lin)/nfft*fs)/pRef));

fprintf('PSD estimated intensity for HF= %.1f dB SPL\n', 20*log10(sqrt((sum(Plow_lin)+sum(Phigh_lin))/nfft*fs)/pRef));
