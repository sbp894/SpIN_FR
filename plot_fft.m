%[P1,f]=plot_fft(X,Fs,PLOT_VAR)
%P1=amplitudes

function [P1,f]=plot_fft(vecin,Fs,x_semilog,title_str)

if nargin==3
    title_str='fft plot';
elseif nargin==2
    x_semilog=1;
    title_str='fft plot';
end

L = length(vecin);
Y = fft(vecin);

P2 = abs(Y/L);
P1 = P2(1:ceil(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
P1=20*log10(P1);
f =linspace(0,Fs/2,length(P1));
if x_semilog
    semilogx(f,P1);
%     xlim([.1 Fs/2]);
else
    plot(f,P1);
end

title(title_str);
xlabel('f (Hz)');
ylabel('20*log10(|P1(f)|), dB');
