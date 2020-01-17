% assuming s_data_env already in workspace
function data_out=temp_remove_lines_chronux(s_data_env, fs_data, plotVar)
clf
s_data_env=s_data_env(:);
fs=fs_data;
params.Fs=fs_data;
params.tapers= [5 9];
N=length(s_data_env );
p=1-1/N;
[~,~,~,~,~,~,params]=getparams(params);


% [~,~,freqs,~,~]=fitlinesc(s_data_env,params,p/N,'n',[]);
nfft=2^nextpow2(N);

[Fval,~,f,~] = ftestc(s_data_env,params,p,'n');

fTest= 75:75:1e3;
freq_out= return_closest_max(f, Fval, fTest, ceil(1/(fs/nfft))); % ceil(1/(fs/nfft)): To check in a 1 Hz range
[datafit,~,~,~,~]=fitlinesc(s_data_env,params,p/N,'n',freq_out);

datan=s_data_env-datafit;

if plotVar
    rmlinesc(s_data_env,params, p, 'n');
    [S2,f]=mtspectrumc(detrend(datan),params);
    subplot(325); hold on;
    plot(f, 10*log10(S2), 'k-','LineWidth', 1.5);
    xlim([10 1e3])
    
    subplot(321); xlim([10 1e3])
    subplot(323); xlim([10 1e3])
end

data_out = datan;