%%
clear;
close all;
clc;

%%
nReps2Avg=200;
nSACs=500;
nSACs2Avg=100;
nPSDs=200;


%%
ChinID=313;
DataRepository='D:\Study Stuff\Matlab\NelData\FFRMATDataResampled\';

checkDIR=dir(sprintf('%s*Q%d*',DataRepository,ChinID));
if isempty(checkDIR)
    error('No such directory for animal number %d',ChinID);
elseif length(checkDIR)~=1
    error('Multiple directories. Change!');
else
    DataDir=[DataRepository checkDIR.name filesep];
end

%%
allfiles=dir([DataDir 'p*FFR*']);
fileVar=4;
load([DataDir allfiles(fileVar).name]);
fprintf('%s\n',allfiles(fileVar).name);
avg_pos=data.AD_Data.AD_Avg_PO_V{1};
avg_neg=data.AD_Data.AD_Avg_NP_V{1};
all_rep=data.AD_Data.AD_All_V;
all_pos=all_rep(1:2:length(all_rep));
all_neg=all_rep(2:2:length(all_rep));

tStart=10e-3;
tEnd=data.Stimuli.slow.duration_ms/1e3;
Fs=data.Stimuli.FsNew;
inds=floor(tStart*Fs:tEnd*Fs);
tFFR=inds/Fs;

%% Parse fname to see what stimulus was played
if data.Stimuli.NoiseType
    StimDir='stimSetFluctuating\';
else
    StimDir='stimSetStationary\';
end
stimfName=[StimDir data.Stimuli.filename];
[dataStim, FsStim]=audioread(stimfName);
tStim=(1:length(dataStim))/FsStim;

%% Analyses using AVGs
avg_pos=avg_pos(inds);
avg_neg=avg_neg(inds);

figure(1);
subplot(4,1,1); plot(tStim, dataStim); ...
    hold on; plot(tStim, -dataStim); ...
    legend('POS', 'NEG'); title('Stimuli');
subplot(4,1,2:4); ...
    hold on; xlabel('Time (sec)'); ...
    plot(tFFR, avg_pos); ...
    plot(tFFR, avg_neg); ...
    legend('POS','NEG'); title('FFR');

%%
avg_sum=avg_pos+avg_neg;
avg_sum=avg_sum-mean(avg_sum);
avg_dif=(avg_pos-avg_neg)/2;
avg_dif=avg_dif-mean(avg_dif);

figure(2);
subplot(5,1,1); plot(tStim, dataStim); ...
    hold on; plot(tStim, -dataStim); ...
    legend('POS', 'NEG'); title('Stimuli');

subplot(5,1,2:3); hold on; xlabel('Time (sec)');
plot(tFFR, avg_sum); plot(tFFR, avg_dif);
legend('SUM','DIFF'); title('FFR');

subplot(5,1,4:5); hold on;
plot_fft(avg_sum, Fs);
plot_fft(avg_dif, Fs);
legend('SUM','DIFF'); title('FFR');


%%
CFs=[250 500 1e3 2e3];
filtOutPos=zeros(length(CFs), length(avg_pos));
filtOutNeg=zeros(length(CFs), length(avg_pos));
hilbertOut=zeros(length(CFs), length(avg_pos));
hilbertAmp=zeros(length(CFs), ceil(length(avg_pos)/2)+1);

sac_pos=zeros(length(CFs), 2*length(avg_pos)-1);
sac_neg=zeros(length(CFs), 2*length(avg_pos)-1);
xpac=zeros(length(CFs), 2*length(avg_pos)-1);
delays=(-length(avg_pos)+1:length(avg_pos)-1)/Fs;

figure(3);
subplot(1+length(CFs),1,1);  plot(tStim, dataStim); ...
    hold on; plot(tStim, -dataStim); ...
    legend('POS', 'NEG'); title('Stimuli');

parfor cfVar=1:length(CFs)
    N=30;
    F3dB1=CFs(cfVar)/sqrt(2);
    F3dB2=CFs(cfVar)*sqrt(2);
    Ast=100;
    D=fdesign.bandpass('N,F3dB1,F3dB2,Ast', N,F3dB1,F3dB2,Ast, Fs);
    dmethods=designmethods(D);
    Hd = design(D,'default');
    %     freqz(Hd,2^12, 'whole', Fs);
    
    filtOutPos(cfVar, :)=filter(Hd, avg_pos);
    filtOutNeg(cfVar, :)=filter(Hd, avg_neg);
    
    hilbertOut(cfVar, :)=(abs(hilbert(filtOutPos(cfVar, :)))+abs(hilbert(filtOutNeg(cfVar, :))))/2;
   
    %%    
    sac_pos(cfVar, :)=xcorr(filtOutPos(cfVar, :));
    sac_neg(cfVar, :)=xcorr(filtOutNeg(cfVar, :));
    xpac(cfVar, :)=xcorr(filtOutPos(cfVar, :), filtOutNeg(cfVar, :));
end

filtOutSum=filtOutPos+filtOutNeg;
filtOutDif=(filtOutPos-filtOutNeg)/2;

sac=(sac_pos+sac_pos)/2;
SUMCOR=sac+xpac;
DIFCOR=(sac-xpac)/2;
figure(4); figure(5);

modFreqs=[1 2 4 8 16 32 64];
MODp=nan(length(CFs), length(modFreqs));

for cfVar=1:length(CFs)
    figure(3);
    subplot(1+length(CFs),1,1+cfVar);  plot(tFFR, filtOutPos(cfVar, :)); hold on; ...
        plot(tFFR, filtOutNeg(cfVar, :)); ylabel(['CF= ' num2str(CFs(cfVar))]); ...
        legend('POS', 'NEG'); title('FFR Data');
    
    figure(4);
    subplot(length(CFs), 2, cfVar*2-1); ...
        plot(delays, SUMCOR(cfVar,:)); ...
        ylabel(['CF= ' num2str(CFs(cfVar))]); ...
        title('SUMCOR');
    subplot(length(CFs), 2, cfVar*2); ...
        plot_fft(SUMCOR(cfVar,:), Fs, 1);
    
    figure(5);
    subplot(length(CFs), 2, cfVar*2-1); ...
        plot(tFFR, hilbertOut(cfVar,:));
    subplot(length(CFs), 2, cfVar*2); ...
        [hilbertAmp(cfVar,:),f]=plot_fft(hilbertOut(cfVar,:), Fs, 1);
    
    PSDcur=hilbertAmp(cfVar,:).^2;
    for modVar=1:length(modFreqs)
       curIndStart= dsearchn(f, modFreqs(modVar)/sqrt(2));
       curIndEnd= dsearchn(f, modFreqs(modVar)*sqrt(2));
       MODp(cfVar, modVar)=sum(PSDcur(curIndStart:curIndEnd));
    end
end


