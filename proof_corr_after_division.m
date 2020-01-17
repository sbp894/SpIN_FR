clear;
clc;

fig_save_dir= sprintf('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/');
if ~isfolder(fig_save_dir)
    mkdir(fig_save_dir);
end

figure(1);
clf;
fSize= 12;


%% Case 1: positive random numbers

x1= 20*rand(500,1);
y1= .7*x1 + randn(size(x1));
x2= 20*rand(500,1);
y2= .9*x2 + randn(size(x1));

subplot(221)
hold on
scatter(x1, y1)
scatter(x2, y2)
title('Case 1: Intercepts=0');
legend('y1 vs x1', 'y2 vs x2', 'Location', 'northwest')
xlabel('X')
ylabel('Y')

subplot(223)
scatter(log(x1./x2), log(y1./y2), 'k')
xlabel('log(x1./x2)');
ylabel('log(y1./y2)');

% Case 2: Random numbers (+ve and -ve)
x1= 20*rand(500,1);
y1= .7*x1 + randn(size(x1)) - 10;
x2= 20*rand(500,1);
y2= .9*x2 + randn(size(x1)) - 10;

subplot(222)
hold on
scatter(x1, y1)
scatter(x2, y2)
title('Case 2: Intercepts !=0');
xlabel('X');
ylabel('Y');

subplot(224)
scatter(log(x1./x2), log(y1./y2), 'k')
xlabel('log(x1./x2)');
ylabel('log(y1./y2)');

set(findall(gcf,'-property','FontSize'),'FontSize',fSize)
set(gcf, 'Units', 'inches', 'Position', [1 1 12 8]);
saveas(gcf, [fig_save_dir 'corr_hypothesis'], 'png');