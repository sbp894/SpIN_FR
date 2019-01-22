close all;
clear;
clc;
clf;

colSpread= .2;
co(1:3,:)= [1-colSpread*rand(3,1) colSpread*rand(3,1) colSpread*rand(3,1)];
co(4:7,:)= [colSpread*rand(4,1) colSpread*rand(4,1) 1-colSpread*rand(4,1)];

lw=2;
lw2= 3;
count=0;

allChins= [191 366 367 365 368 369 370];
allSPL= [55 70 82];
CodesDir= '/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/';
DataDir= '/media/parida/DATAPART1/Matlab/ExpData/MatData/';
OutDir= 'Figure_Out/';

bxx= nan(length(allSPL), 1);
cxx= nan(length(allSPL), 1);

for splVar= 1:3
    figure(splVar);
    clf;
end

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
        
        
        plot_data_pos_neg=0;
        if plot_data_pos_neg
            figure
            plot(t_data, data_pos_bp)
            hold on
            plot(t_data, data_neg_bp)
            A= .1;
            plot(t_data, A+ data_env);
            plot(t_data, -A+ data_tfs);
            title(ttlStr)
        end
        
        
        
        N_bp_half= 10;
        HalfPowerFrequency1=90;
        HalfPowerFrequency2=130;
        curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
            'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
            'SampleRate',fs_data);
        
        env_f0= filter(curFilt, data_env);
        tfs_f0= filter(curFilt, data_tfs);
        
        ttlStr= sprintf('Q%d, SPL=%.0f, %d-%d Hz band', chinID, allSPL(splVar), HalfPowerFrequency1, HalfPowerFrequency2);
        
        plot_data_env_tfs= 0;
        if plot_data_env_tfs
            figure
            plot(t_data, env_f0)
            hold on
            plot(t_data, tfs_f0)
            title(ttlStr)
        end
        
        segRes= 100e-3;
        fracMove= .5;
        env_dBSPL= gen_get_spl_vals(env_f0, fs_data, segRes, fracMove);
        [tfs_dBSPL, timeVals]= gen_get_spl_vals(tfs_f0, fs_data, segRes, fracMove);
        
        
        plot_rms= 1;
        if plot_rms
            figure(splVar);
            %             yyaxis left
            subplot(211);
            bxx(splVar)= gca;
            hold on
            plot(timeVals, env_dBSPL, '-d', 'color', co(chinVar, :), 'linew', lw);
            plot(timeVals, tfs_dBSPL, '--o', 'color', co(chinVar, :), 'linew', lw);
            title(ttlStr);
            grid on
            
            %             yyaxis right
            subplot(212);
            hold on;
            cxx(splVar)= gca;
            plot(timeVals, env_dBSPL-tfs_dBSPL, '-', 'color', co(chinVar, :), 'linew', lw2);
            grid on
            xlabel('time (sec)');
            
            %             ylim([30 52]);
        end
        
        % figure; scatter(env_rms(1:end-2), tfs_rms(1:end-2))
        
        plot_psd= 0;
        if plot_psd
            count= count+1;
            figure(3+count)
            NW= 2;
            clf;
            [Pxx_env_dB, ~, px_env]= plot_dpss_psd(data_env, fs_data, 'NW', NW, 'plot', true); %#ok<*ASGLU>
            hold on
            [Pxx_tfs_dB, ~, px_env]= plot_dpss_psd(data_tfs, fs_data, 'NW', NW, 'plot', true); %#ok<*ASGLU>
            title(ttlStr)
            axx(count) = gca;
            set(gca, 'xscale', 'lin');
            xlim([10 3e3]);
        end
    end
end
cd (CodesDir);
if exist('axx', 'var')
    linkaxes(axx);
end
if ~prod(isnan(bxx))
    linkaxes(bxx);
end
if ~prod(isnan(cxx))
    linkaxes(cxx);
end
%%
for splVar= 1:3
    figure(splVar);
    lg= nan(3*length(allChins), 1);
    lg_short= nan(6, 1);
    lgStr= cell(3*length(allChins), 1);
    for chinVar= 1:length(allChins)
        chinID= allChins(chinVar);
        lg(3*chinVar-2)= plot(nan, nan, '-', 'color', co(chinVar,:), 'linew', lw);
        lg(3*chinVar-1)= plot(nan, nan, '--', 'color', co(chinVar,:), 'linew', lw);
        lg(3*chinVar)= plot(nan, nan, '-', 'color', co(chinVar,:), 'linew', lw2);
        lgStr{3*chinVar-2}= [ num2str(chinID) '-ENV'];
        lgStr{3*chinVar-1}= [ num2str(chinID) '-TFS'];
        lgStr{3*chinVar}= [ num2str(chinID) '-Ratio'];
    end
    
    figure(splVar);
    subplot(211);
    lg_short(1)= plot(nan, nan, '-', 'color', get(lg(end), 'color'), 'linew', lw);
    lg_short(2)= plot(nan, nan, '--', 'color', get(lg(end), 'color'), 'linew', lw);
    lg_short(3)= plot(nan, nan, '-', 'color', get(lg(1), 'color'), 'linew', lw);
    lg_short(4)= plot(nan, nan, '--', 'color', get(lg(1), 'color'), 'linew', lw);
    legend(lg_short(1:4), {'NH-ENV', 'NH-TFS', 'HI-ENV', 'HI-ENV'}, 'location', 'southwest');
    ylabel('Power in EFR');
    
    subplot(212);
    lg_short(5)= plot(nan, nan, '-', 'color', get(lg(end), 'color'), 'linew', lw2);
    lg_short(6)= plot(nan, nan, '-', 'color', get(lg(1), 'color'), 'linew', lw2);
    legend(lg_short(5:6), {'NH-DIFF', 'HI-ENV'}, 'location', 'southwest');
    ylabel('ENV-TFS power (diff)');

    %     figure(splVar);
    %     subplot(211);
    %     inds_1= setxor(1:3*length(allChins), 3:3:3*length(allChins));
    %     legend(lg(inds_1), lgStr(inds_1), 'location', 'southwest');
    %     figure(splVar);
    %     subplot(212);
    %     inds_2= 3:3:3*length(allChins);
    %     legend(lg(inds_2), lgStr(inds_2), 'location', 'southwest');
    
    set(gcf, 'units', 'normalized', 'position', [.1 .1 .8 .8]);
    figName= sprintf('%sspl%d_f0_%dto%d_Hz', OutDir, allSPL(splVar), HalfPowerFrequency1, HalfPowerFrequency2);
    saveas(gcf, figName, 'tiff');
end