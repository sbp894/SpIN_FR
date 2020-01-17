function latency_val= check_acosutic_delay(fName)
[sig, fs_sig]= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_06_01-Q346_AN_NH/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');

data= load(fName);
data= data.data;

fs_data= round(data.Stimuli.RPsamprate_Hz);
data_pos= data.AD_Data.AD_Avg_PO_V{1}';
data_neg= data.AD_Data.AD_Avg_NP_V{1}';

fs= fs_sig;
data_pos= gen_resample(data_pos, fs_data, fs);
data_neg= gen_resample(data_neg, fs_data, fs);

sig= [sig; zeros(length(data_pos)-length(sig), 1)];

t_data=  (1:length(data_pos))/fs;
t_sig= (1:length(sig))/fs;

[ccf_pos, delay_pos]= xcorr(data_pos, sig, 'coeff');
[ccf_neg, delay_neg]= xcorr(data_neg, -sig, 'coeff');

[~, ind_max_pos]= max(ccf_pos);
[~, ind_max_neg]= max(ccf_neg);

%%
figure(1);
clf;
ax(1)= subplot(221);
yyaxis left;
plot(t_sig*1e3, sig);
grid on;

yyaxis right;
plot(t_data*1e3, data_pos);
xlabel('t in ms');

ax(2)= subplot(222);
yyaxis left;
plot(t_sig*1e3, -sig);
xlabel('t in ms');
grid on;

yyaxis right;
plot(t_data*1e3, data_neg);

linkaxes(ax, 'x');

bx(1)= subplot(223);
hold on;
grid on;
plot(delay_pos/fs*1e3, ccf_pos);
plot(delay_pos(ind_max_pos)/fs*1e3, ccf_pos(ind_max_pos), '*');
xlabel('t in ms');

bx(2)= subplot(224);
hold on;
grid on;
plot(delay_neg/fs*1e3, ccf_neg);
plot(delay_pos(ind_max_neg)/fs*1e3, ccf_pos(ind_max_neg), '*');
xlabel('t in ms');

linkaxes(bx, 'x');
xlim([-10 10]);

latency_val.pos= delay_pos(ind_max_pos)/fs;
latency_val.neg= delay_neg(ind_max_neg)/fs;