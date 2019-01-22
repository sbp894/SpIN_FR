clear;
clc;

fs= 100e3;
fStart= 500;
fEnd= 10e3;
waitDur= 1.5;
dur= 1;
amp=1;
rampTime= 25e-3;

rng(0);

f0_1= 400;
f0_2= 200;

freq1= fStart:f0_1:fEnd;
amp1= hamming(length(freq1)); 
x1= gen_ramp_signal(create_sinusoid(freq1, fs, dur, amp1), fs, rampTime);

freq2= f0_1/2+fStart:f0_1:fEnd;
amp2= hamming(length(freq2));
x2= gen_ramp_signal(create_sinusoid(freq2, fs, dur, amp2), fs, rampTime);

freq3= fStart:f0_2:fEnd;
amp3= hamming(length(freq3));
x3= gen_ramp_signal(create_sinusoid(freq3, fs, dur, amp3), fs, rampTime);

figure(1);
clf;
subplot(211);
hold on;
plot_dpss_psd(x1, fs);
plot_dpss_psd(x2, fs);
subplot(212);
plot_dpss_psd(x3, fs);

t=(1:length(x1))/fs;

sound(x1, fs);
pause(waitDur*dur);
sound(x2, fs);
pause(waitDur*dur);
sound(x3, fs);