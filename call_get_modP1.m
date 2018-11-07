% clear;
close all;
clc;


%%
lw=2;
cols='brkm';

ChinID=313;


fileVar=2;
[MODp_SN, ~, ~]=get_modP(ChinID, fileVar);

fileVar=1;
[MODp_N, modFreqs, CFs]=get_modP(ChinID, fileVar);

SNRenvAll=(MODp_SN-MODp_N)./MODp_N;
SNRenvAll(SNRenvAll<0)=0;

close all; 
figure(1); hold on;
figure(2); hold on;

for cfVar=1:size(MODp_SN,1)
   figure(1); ...
       plot(MODp_SN(cfVar,:), [cols(cfVar) '-o'], 'linewidth',lw); ...
       plot(MODp_N(cfVar,:), [cols(cfVar) '--d'], 'linewidth',lw); 
   figure(2); ...
       plot(SNRenvAll(cfVar,:), [cols(cfVar) '-o'], 'linewidth',lw); ...   
end

figure(1); 
set(gca, 'fontsize', 16, 'fontweight', 'bold', 'xticklabel',  strread(num2str(modFreqs),'%s'));
xlabel('mod freq (hz)'); legend(strread(num2str(CFs),'%s'));

figure(2); 
set(gca, 'fontsize', 16, 'fontweight', 'bold', 'xticklabel',  strread(num2str(modFreqs),'%s'));
xlabel('mod freq (hz)');  legend(strread(num2str(CFs),'%s'));


SNRenv2=20*log10(sum(sum(SNRenvAll.^2)))