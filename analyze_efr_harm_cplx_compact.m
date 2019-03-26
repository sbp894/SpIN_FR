clear;
clc;
clf;

colSpread= .2;

lw=2;
lw2= 3;
mrkSize= 16;
count=0;
fSize= 16;

allChins= [191 366 367 365 368 369 370];
co(ismember(allChins, [191 366 367]),:)= [1-colSpread*rand(3,1) colSpread*rand(3,1) colSpread*rand(3,1)];
co(ismember(allChins, [365 368 369 370]),:)= [colSpread*rand(4,1) colSpread*rand(4,1) 1-colSpread*rand(4,1)];

allSPL= [55 70 82];
CodesDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/';
DataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
outFigDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/Figure_Out/EFR_hrm_cmplx/';


[hrm_cplx_sig, fs_sig]= audioread('/media/parida/DATAPART1/Matlab/Design_Exps_NEL/create_harmonic_complex_ffr_hf_vs_lf/LF_HF_CMPLX_HRMNC_Stimuli/LFlow2high_complex.wav');
bp_hf_audio_Filt= get_filter(fs_sig, 500, 3.1e3);
bp_lf_audio_Filt= get_filter(fs_sig, 80, 500);
hrm_cplx_sig_hf= filtfilt(bp_hf_audio_Filt, hrm_cplx_sig);
hrm_cplx_sig_lf= filtfilt(bp_lf_audio_Filt, hrm_cplx_sig);

nSProws= 2;
nSPcols= length(allSPL);
axx= nan(nSPcols, 1);
bxx= nan(nSPcols, 1);


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
        
        fs_data= x.Stimuli.RPsamprate_Hz;
        tMax= x.Stimuli.FFR_Gating.duration_ms/1e3;
        dataLen= round(fs_data*tMax);
        t_data= (1:dataLen)/fs_data;
        
        
        data_pos= data_pos(1:dataLen);
        data_neg= data_neg(1:dataLen);
        t_data= t_data(1:dataLen);
        
        
        N_bp_half= 4;
        HalfPowerFrequency1=50;
        HalfPowerFrequency2=3.2e3;
        curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
            'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
            'SampleRate',x.Stimuli.RPsamprate_Hz);
        
        data_pos_bp= filter(curFilt, data_pos);
        data_neg_bp= filter(curFilt, data_neg);
        data_env= (data_pos_bp + data_neg_bp)/2;
        data_tfs= (data_pos_bp - data_neg_bp)/2;
        
        N_bp_half= 10;
        HalfPowerFrequency1=90;
        HalfPowerFrequency2=460;
        curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
            'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
            'SampleRate',fs_data);
        
        env_f0= filter(curFilt, data_env);
        tfs_f0= filter(curFilt, data_tfs);
        
        ttlStr= sprintf('%.0f dB SPL', allSPL(splVar));
        
        segRes= 200e-3;
        fracMove= .5;
        env_dBSPL= gen_get_spl_vals(env_f0, fs_data, segRes, fracMove);
        [tfs_dBSPL, timeVals]= gen_get_spl_vals(tfs_f0, fs_data, segRes, fracMove);
        sig_hf_spl= gen_get_spl_vals(hrm_cplx_sig_hf, fs_sig, segRes, fracMove);
        sig_lf_spl= gen_get_spl_vals(hrm_cplx_sig_lf, fs_sig, segRes, fracMove);
        
        plot_rms= 1;
        if plot_rms
            axx(splVar)=subplot(nSProws, nSPcols, splVar);
            hold on
            plot(timeVals, env_dBSPL, '-d', 'color', co(chinVar, :), 'linew', lw, 'markersize', mrkSize);
            grid on;
            title(ttlStr);
            set(gca, 'fontsize', fSize);
            if splVar~=1
                set(gca, 'yticklabel', '');
            end
            
            bxx(splVar)=subplot(nSProws, nSPcols, nSPcols+splVar);
            hold on
            plot(timeVals, tfs_dBSPL, '-o', 'color', co(chinVar, :), 'linew', lw, 'markersize', mrkSize);
            grid on
            set(gca, 'fontsize', fSize);
            if splVar~=1
                set(gca, 'yticklabel', '');
            end
        end
    end
end

subplot(nSProws, nSPcols, nSProws*nSPcols-1);
xlabel('Time (sec)');
subplot(nSProws, nSPcols, 1);
ylabel('$ENV_{Power}$', 'interpreter', 'latex');
subplot(nSProws, nSPcols, (nSProws-1)*nSPcols+1);
ylabel('$TFS_{Power}$', 'interpreter', 'latex');
subplot(nSProws, nSPcols, nSProws*nSPcols);
dummy_lines(1)= plot(nan, nan, 'b', 'linew', lw);
dummy_lines(2)= plot(nan, nan, 'r', 'linew', lw);
legend(dummy_lines, 'NH', 'HI', 'fontsize', 10, 'location', 'northwest');


cd (CodesDir);
if ~prod(isnan(axx))
    linkaxes(axx);
end
if ~prod(isnan(bxx))
    linkaxes(bxx);
end

set(gcf, 'units', 'inches', 'position', [1 1 11 6]);
fName= 'hrm_cplx_all_spl';
saveas(gcf, [outFigDir fName], 'png');

function bpFilt= get_filter(fs_data, HalfPowerFrequency1, HalfPowerFrequency2)
N_bp_half= 4;
% HalfPowerFrequency1=500;
% HalfPowerFrequency2=3.1e3;

bpFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs_data);
end
