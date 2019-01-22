% clear;
clc;
clf


f0_base= 100;
dur= .5;
f0_test_init= 105 + 4*rand(1);
f0_steps= [2 .5 .2];
nReversals= [2 3 5];

fs= 50e3;

low_freq.nHarmStart= 3;
low_freq.nHarmEnd= 5;
low_freq.base_freqs= f0_base*(low_freq.nHarmStart:low_freq.nHarmEnd);
low_freq.base_sig= create_harmCplx(f0_base, fs, low_freq.nHarmStart:low_freq.nHarmEnd, dur);

stepData= struct('freqs', nan, 'f0_test', nan, 'f0_step', nan, 'resp', nan, 'trueInterval', nan, 'UserInInterval', nan);
count=1;


sigs2play= cell(3,1);
sigs2play{1}= low_freq.base_sig;

stopFlag= 0;
cur_f0_Step= f0_steps(1);
f0_test= f0_test_init;
while ~stopFlag
    curFreqs= f0_test*(low_freq.nHarmStart:low_freq.nHarmEnd);
    sigs2play{2}= create_harmCplx(f0_base, fs, low_freq.nHarmStart:low_freq.nHarmEnd, dur);
    sigs2play{3}= create_harmCplx(f0_test, fs, low_freq.nHarmStart:low_freq.nHarmEnd, dur);
    play_order= randsample(3, 3);
    for playVar= 1:length(play_order)
%         sound(sigs2play{play_order(playVar)}, fs);
        %         pause(1.5*dur);
    end
    
    %     usr_ans= str2double(questdlg('Which one diff?', 'Select', '1', '2', '3', '2'));
    errorVal= 1./(1.75*cur_f0_Step*randn(1))
    usr_ans= dsearchn((1:3)', find(play_order==3)+ errorVal);
    %     pause(dur);
    
    stepData(count).freqs= curFreqs;
    stepData(count).f0_test= f0_test;
    stepData(count).f0_step= cur_f0_Step;
    stepData(count).trueInterval= find(play_order==3);
    stepData(count).UserInInterval= usr_ans;
    
    
    curStepInds= [stepData.f0_step]== cur_f0_Step;
    if usr_ans == find(play_order==3)
        fprintf('Bingo\n');
        % Correct
        stepData(count).resp= true;
        ttlStr= 'correct';
        
        % Check how many correct cons.
        if count>2
            if all([stepData(count-2:count).resp]) && length(unique([stepData(count-2:count).f0_test]))==1
                f0_test= f0_test - cur_f0_Step;
                if f0_test<f0_base
                    f0_test=f0_base;
                end
            end
        end
        
    else
        % Wrong
        stepData(count).resp= false;
        ttlStr= 'wrong';
        fprintf('Boo\n');
        f0_test= f0_test+ cur_f0_Step;
    end
    
    
    
    % Check the number of reversals for this step
    check_f0_test_data= [[stepData(curStepInds).f0_test] f0_test];
    if length(check_f0_test_data)>2
        nRev_test= numel(findpeaks(check_f0_test_data)) + numel(findpeaks(-check_f0_test_data));
        nRev_init= nReversals(cur_f0_Step==f0_steps);
        if nRev_test==nRev_init
            % Reversals done. Either change step size or stop
            if cur_f0_Step~=f0_steps(end)
                % Change step size
                cur_f0_Step= f0_steps(find(f0_steps==cur_f0_Step)+1);
            else
                % Stop
                break;
            end
        end
    end
    
    figure(1);
    plot([[stepData.f0_test] f0_test], 'd-');
    title(ttlStr);
    
    count= count+1;
end

for stepVar= 1:length(f0_steps)
    cur_f0_Step= f0_steps(stepVar);
    
    
    
    stepData(count).f0_test= cur_f0_Step;
    stepData(count).f0_step= cur_f0_Step;
end


high_freq.nHarmStart= 10*low_freq.nHarmStart;
high_freq.nHarmEnd= 10*low_freq.nHarmEnd;
high_freq.base_freqs= f0_base*(high_freq.nHarmStart:high_freq.nHarmEnd);



function sig= create_harmCplx(f0, fs, nHarm, dur)
ramp= .05;
snr= 20;
sig= create_sinusoid(f0*nHarm, fs, dur);
sig= gen_ramp_signal(sig, fs, ramp);
sig= create_noisy_signal(sig, snr, 'pink');
end