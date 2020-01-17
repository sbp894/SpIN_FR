function out_data=remove_artifact_ffr(in_data, fs_data, plotVar, fMax, x_output_spread)

if ~exist('x_output_spread', 'var')
    x_output_spread= 0;
end

if ~exist('plotVar', 'var')
    plotVar= 0;
end

if ~exist('fMax', 'var')
    fMax= 500;
end

if size(in_data,1)==1 % means row vector
    rowV1_colV0= 1;
else
    rowV1_colV0= 0;
end

in_data=in_data(:);
fs=fs_data;
params.Fs=fs_data;
params.pad=2;
params.tapers= [5 9];
N=length(in_data );
p=1-1/N;
[~,~,~,~,~,~,params]=getparams(params);


% [~,~,freqs,~,~]=fitlinesc(s_data_env,params,p/N,'n',[]);
nfft=2^nextpow2(N);

[Fval,~,f,~] = ftestc(in_data,params,p,'n');

fTest= sort([75:75/2:fMax 180]);
freq_out= return_closest_max(f, Fval, fTest, ceil(2/(fs/nfft)), x_output_spread); % ceil(1/(fs/nfft)): To check in a 1 Hz range
[datafit,~,~,~,~]=fitlinesc(in_data,params,p/N,'n',freq_out);

out_data=in_data-datafit;
xtick_vals= [50 100 500];
xtick_labs= cellfun(@(x) num2str(x), num2cell(xtick_vals), 'UniformOutput', false);
if plotVar
    hold on;
    %     yrange=32;
    plot_dpss_psd(in_data, fs_data, 'nw', params.tapers(1));
    plot_dpss_psd(out_data, fs_data, 'nw', params.tapers(1));
    xlim([25 fMax+50]);
    set(gca, 'XTick', xtick_vals, 'XTickLabel', xtick_labs);
end

if rowV1_colV0
    out_data= out_data';
end