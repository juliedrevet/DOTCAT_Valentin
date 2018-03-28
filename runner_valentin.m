%%
close all
clear all java
close all hidden
clc


subj = 999;

% customizable options:
runbgncalib = true;
runmidcalib = true;
runendcalib = true;
start_blck  = 5;

%% create data folder for subject
foldname = sprintf('./Data/S%02d',subj);
if ~exist(foldname,'dir')
    mkdir(foldname);
end

diary(sprintf('./Data/S%02d/S%02d_Log_%s.txt',subj,subj,datestr(now,'yyyymmdd-HHMM')));

%% RUN BEGINNING CALIBRATION
fprintf('*** BEGINNING CALIBRATION ***\n\n')
fprintf('load calibration file...')
calibration = importdata('./Calibration/DOTCAT_calibration_raw_valentin.mat');
fprintf('done!\n\n')

fprintf('Press any key to continue.\n')
pause
[calibration1,~,~] = run_calibration(subj,calibration,true); % with mini training!
fprintf('\t min / max proportion of Color 1: \t%d / %d\n\n',round(calibration1.rslt.range(1)*100),round(100*calibration1.rslt.range(2)))


%% RUN EXPERIMENT 1st HALF WITHOUT PRACTICE
fprintf('*** EXPERIMENT 1st Half ***\n\n')
expe  = gen_expe(subj);
shape = shuffle_shape(expe); % create subject shape combinations
expe_raw       = [];
expe_raw.expe  = expe;
expe_raw.shape = shape;
% save raw expe structure in case of premature termination
fname = sprintf('DOTCAT_S%02d_%s_raw',subj,datestr(now,'yyyymmdd-HHMM'));
fname = fullfile(foldname,fname);
save([fname,'.mat'],'expe_raw');


fprintf('Press any key to continue.\n')
pause
[expe,~] = run_expe(subj,calibration1,expe_raw,start_blck,8);
start_blck = 9;


%% RUN MID-CALIBRATION
fprintf('*** MID-CALIBRATION ***\n\n')
fprintf('Press any key to continue.\n')
pause
[calibration2,~,~] = run_calibration(subj,calibration,false); % no training
fprintf('\t min / max proportion of Color 1: \t%d / %d\n\n',round(calibration2.rslt.range(1)*100),round(100*calibration2.rslt.range(2)))

%% RUN EXPERIMENT 2nd HALF
fprintf('*** EXPERIMENT 2nd Half ***\n\n')
fprintf('Press any key to continue.\n')
pause
[expe,~] = run_expe(subj,calibration2,expe_raw,start_blck,12);

%% RUN END CALIBRATION
fprintf('*** END CALIBRATION ***\n\n')
fprintf('Press any key to continue.\n')
pause
[calibration3,aborted,errmsg] = run_calibration(subj,calibration,false); % no training


%% PLOT CALIB
plot_calib(calibration1);
suptitle(sprintf('Subject %d : beginning calibration',subj))
savefig(sprintf('./Data/S%02d/calibration1',subj))
plot_calib(calibration2);
suptitle(sprintf('Subject %d : mid-calibration',subj))
savefig(sprintf('./Data/S%02d/calibration2',subj))
plot_calib(calibration3);
suptitle(sprintf('Subject %d : end calibration',subj))
savefig(sprintf('./Data/S%02d/calibration3',subj))
%
plot_calib(calibration1, false,true);
plot_calib(calibration2,true,true);
plot_calib(calibration3,true,true);
savefig(sprintf('./Data/S%02d/calibration_all',subj))

fprintf('\nBeginning Calibration:\n')
fprintf('\t min / max proportion of Color 1: \t%d / %d\n\n',round(calibration1.rslt.range(1)*100),round(100*calibration1.rslt.range(2)))
fprintf('Mid Calibration:\n')
fprintf('\t  min / max proportion of Color 1: \t%d / %d\n\n',round(100*calibration2.rslt.range(1)),round(100*calibration2.rslt.range(2)))
fprintf('End Calibration:\n')
fprintf('\t  min / max proportion of Color 1: \t%d / %d\n\n',round(100*calibration3.rslt.range(1)),round(100*calibration3.rslt.range(2)))



diary('off')