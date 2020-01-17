clear;
clc;
clf;

CodesDirs= {'/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/helper', ...
    '/media/parida/DATAPART1/Matlab/MWcentral/chronux_2_11/spectral_analysis/continuous'};
addpath(CodesDirs{:});


colSpread= .2;
nhChins= [365 368 369 370];
hiChins= [366 367];
allChins= [nhChins hiChins]; % 191
co(ismember(allChins, hiChins),:)= [1-colSpread*rand(numel(hiChins),1) colSpread*rand(numel(hiChins),1) colSpread*rand(numel(hiChins),1)];
co(ismember(allChins, nhChins),:)= [colSpread*rand(numel(nhChins),1) colSpread*rand(numel(nhChins),1) 1-colSpread*rand(numel(nhChins),1)];

lw=2;
lw2= 3;
mrkSize= 6;
count=0;
fSize= 20;

allSPL= [55 70 82];
CodesDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/';
DataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
OutFigDir= [CodesDir 'Figure_Out/EFR_hrm_cmplx/IndvPlots/'];

if ~isfolder(OutFigDir)
    mkdir(OutFigDir);
end

bxx= nan(length(allSPL), 1);
cxx= nan(length(allSPL), 1);

% for splVar= 1:3
%     figure(splVar);
%     clf;
% end

[hrm_cplx_sig, fs_sig]= audioread('/media/parida/DATAPART1/Matlab/Design_Exps_NEL/create_harmonic_complex_ffr_hf_vs_lf/LF_HF_CMPLX_HRMNC_Stimuli/LFlow2high_complex.wav');
bp_hf_audio_Filt= get_filter(fs_sig, 500, 3.1e3);
bp_lf_audio_Filt= get_filter(fs_sig, 80, 500);
hrm_cplx_sig_hf= filtfilt(bp_hf_audio_Filt, hrm_cplx_sig);
hrm_cplx_sig_lf= filtfilt(bp_lf_audio_Filt, hrm_cplx_sig);

for splVar= 1:length(allSPL)
    for chinVar= 1:length(allChins)
        chinID= allChins(chinVar);
        
        
        allDirs= dir([DataDir '*' num2str(chinID) '*EFR*']);
        cd([DataDir allDirs.name]);
        
        allFiles= dir('*EFR*');
        attns= cell2mat(cellfun(@(x) str2double(x(28:end-2)), {allFiles.name}', 'uniformoutput', false));
        [~, sort_ind]= sort(attns, 'descend');
        
        mfilename= allFiles(sort_ind==splVar).name;
        x= load(mfilename);
        x= x.data;
        
        data_pos= x.AD_Data.AD_Avg_PO_V{1};
        data_neg= x.AD_Data.AD_Avg_NP_V{1};
        
        fs_data= 10e3;
        data_pos= gen_resample(data_pos, x.Stimuli.RPsamprate_Hz, fs_data);
        data_neg= gen_resample(data_neg, x.Stimuli.RPsamprate_Hz, fs_data);
        
        plotVar= false;
        data_pos=remove_artifact_ffr(data_pos, fs_data, plotVar);
        data_neg=remove_artifact_ffr(data_neg, fs_data, plotVar);
        
        
        tMax= x.Stimuli.FFR_Gating.duration_ms/1e3;
        dataLen= round(fs_data*tMax);
        t_data= (1:dataLen)/fs_data;
        
        
        data_pos= data_pos(1:dataLen);
        data_neg= data_neg(1:dataLen);
        t_data= t_data(1:dataLen);
        
        
        N_bp_half= 4;
        HalfPowerFrequency1=100;
        HalfPowerFrequency2=3.2e3;
        curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
            'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
            'SampleRate',x.Stimuli.RPsamprate_Hz);
        
        data_pos_bp= filter(curFilt, data_pos);
        data_neg_bp= filter(curFilt, data_neg);
        data_env= (data_pos_bp + data_neg_bp)/2;
        data_tfs= (data_pos_bp - data_neg_bp)/2;
        
        
        if ismember(chinID, nhChins)
            hearing_stat= 'NH';
        else 
            hearing_stat= 'HI';
        end
        ttlStr= sprintf('Q%d | %s |Intensity=%.0f dB SPL, %d-%d Hz band', chinID, hearing_stat, allSPL(splVar), HalfPowerFrequency1, HalfPowerFrequency2);
        
        plot_data_env_tfs= 1;
        if plot_data_env_tfs
            figure(1);
            clf;
            subplot(2,2,1:2);
            plot(t_data, data_env)
            hold on
            plot(t_data, data_tfs)
            title(ttlStr)
        end
        
        
        
        limitFreq= 0;
        if limitFreq
            N_bp_half= 10;
            HalfPowerFrequency1=90;
            HalfPowerFrequency2=460;
            curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
                'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
                'SampleRate',fs_data);
            
            env_f0= filter(curFilt, data_env);
            tfs_f0= filter(curFilt, data_tfs);
        else
            env_f0= data_env;
            tfs_f0= data_tfs;
        end
        
        ttlStr= sprintf('Intensity=%.0f dB SPL, %d-%d Hz band', allSPL(splVar), HalfPowerFrequency1, HalfPowerFrequency2);
        plot_data_env_tfs= 0;
        if plot_data_env_tfs
            figure(1);
            subplot(2,2,1:2);
            plot(t_data, env_f0)
            hold on
            plot(t_data, tfs_f0)
            title(ttlStr)
        end
        set(gca, 'FontSize', fSize);
        
        segRes= 80e-3;
        fracMove= .5;
        env_dBSPL= gen_get_spl_vals(env_f0, fs_data, segRes, fracMove)+20*log10(20e-6);
        [tfs_dBSPL, timeVals]= gen_get_spl_vals(tfs_f0, fs_data, segRes, fracMove);
        tfs_dBSPL= tfs_dBSPL+20*log10(20e-6);
        sig_hf_spl= gen_get_spl_vals(hrm_cplx_sig_hf, fs_sig, segRes, fracMove);
        sig_lf_spl= gen_get_spl_vals(hrm_cplx_sig_lf, fs_sig, segRes, fracMove);
        
        plot_rms= 1;
        if plot_rms
            figure(1);
            subplot(223);
            hold on;
            plot(timeVals, env_dBSPL, '-d', 'color', 'k', 'linew', lw, 'markersize', mrkSize);
            plot(timeVals, tfs_dBSPL, '--o', 'color', 'k', 'linew', lw, 'markersize', mrkSize);
            grid on
            set(gca, 'FontSize', fSize);
            legend('ENV', 'TFS');
            ylabel('short term power (dB)')
        end
        
        
        plot_psd= 1;
        if plot_psd
            figure(1);
            subplot(224)
            xtick_vals= [50 100 300 1e3 3e3];
            xtick_labs= cellfun(@(x) num2str(x), num2cell(xtick_vals), 'UniformOutput', false);
            
            NW= 10;
            [Pxx_env_dB, ~, px_env]= plot_dpss_psd(data_env, fs_data, 'NW', NW, 'plot', true, 'yrange', 75); %#ok<*ASGLU>
            hold on
            [Pxx_tfs_dB, freq, px_tfs]= plot_dpss_psd(data_tfs, fs_data, 'NW', NW, 'plot', false, 'yrange', 75);
            set(px_env, 'linew', 3)
            plot(freq+5, Pxx_tfs_dB, '-', 'linew', 2);
            %             title(ttlStr)
            xlim([50 5e3]);
            ylim([-120 -40])
            set(gca,'XTick', xtick_vals, 'XTickLabel', xtick_labs, 'FontSize', fSize);
            legend('ENV', 'TFS');
        end
        set(gcf, 'Units', 'inches', 'Position', [1 1 16 9]);
        saveas(gcf, [OutFigDir sprintf('Q%d_spl%d', chinID, allSPL(splVar))], 'png');
    end
end
cd (CodesDir);
rmpath(CodesDirs{:});

function bpFilt= get_filter(fs_data, HalfPowerFrequency1, HalfPowerFrequency2)
N_bp_half= 4;
% HalfPowerFrequency1=500;
% HalfPowerFrequency2=3.1e3;

bpFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end
