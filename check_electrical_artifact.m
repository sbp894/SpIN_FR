clear;
clc;



nh.Dir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2019_03_19-Q371_SFR_pink_NH/';
nh.allfiles= dir([nh.Dir '*.mat']);

nh.art_file= nh.allfiles(contains({nh.allfiles.name}, '_artifact')).name;
nh.reg_file= nh.allfiles(contains({nh.allfiles.name}, 'snr_120') & ~contains({nh.allfiles.name}, '_artifact')).name;

temp= load([nh.Dir nh.art_file]);
nh.art_data.pos= temp.data.AD_Data.AD_Avg_PO_V{1};
nh.art_data.neg= temp.data.AD_Data.AD_Avg_NP_V{1};

temp= load([nh.Dir nh.reg_file]);
nh.reg_data.pos= temp.data.AD_Data.AD_Avg_PO_V{1};
nh.reg_data.neg= temp.data.AD_Data.AD_Avg_NP_V{1};

fs_data= temp.data.Stimuli.RPsamprate_Hz;
t_data= (1:length(nh.reg_data.pos))/fs_data;

HalfPowerFrequency1= 70;
HalfPowerFrequency2= 3e3;
N_bp_half= 4;

bpFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);

%
nh.reg_data.pos= filtfilt(bpFilt, nh.reg_data.pos);
nh.reg_data.neg= filtfilt(bpFilt, nh.reg_data.neg);
nh.reg_data.env= (nh.reg_data.pos+nh.reg_data.neg)/2;
nh.reg_data.tfs= (nh.reg_data.pos-nh.reg_data.neg)/2;

%
nh.art_data.pos= filtfilt(bpFilt, nh.art_data.pos);
nh.art_data.neg= filtfilt(bpFilt, nh.art_data.neg);
nh.art_data.env= (nh.art_data.pos+nh.art_data.neg)/2;
nh.art_data.tfs= (nh.art_data.pos-nh.art_data.neg)/2;


%%
figure(1);
clf;

subplot(221);
hold on;
grid on;
plot(t_data, nh.reg_data.pos);
plot(t_data, nh.reg_data.neg);

subplot(222);
hold on;
grid on;
plot(t_data, nh.reg_data.env);
plot(t_data, nh.reg_data.tfs);


yRange= 50;
nfft= round(fs_data);
subplot(223);
hold on;
plot_dpss_psd(nh.reg_data.pos, fs_data, 'nfft', nfft);
[~,~,l2]= plot_dpss_psd(nh.art_data.pos, fs_data, 'nfft', nfft);
set(l2, 'color', 'g', 'LineStyle', '--');
axis tight;
yl= ylim;
ylim([max(yl)-yRange max(yl)]);

subplot(224);
hold on;
plot_dpss_psd(nh.reg_data.env, fs_data, 'nfft', nfft);
plot_dpss_psd(nh.reg_data.tfs, fs_data, 'nfft', nfft);
[~,~,l3]= plot_dpss_psd(nh.art_data.env, fs_data, 'nfft', nfft);
set(l3, 'color', 'g', 'LineStyle', '--');
[~,~,l4]= plot_dpss_psd(nh.art_data.tfs, fs_data, 'nfft', nfft);
set(l4, 'color', 'c', 'LineStyle', '--');
axis tight;
yl= ylim;
ylim([max(yl)-yRange max(yl)]);