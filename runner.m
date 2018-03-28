%%
close all
clear all java
close all hidden
clc

addpath ./Analysis
% ask for subject number and raw experiment structure

prompt={'Subject number (two-digit)','load raw experiment structure?'};
argindlg = inputdlg(prompt,'DOTCAT',1,{'','no'});


if isempty(argindlg)
    error('experiment cancelled!');
end
if isempty(argindlg{1})
    error('experiment cancelled!');
end

subj = str2num(argindlg{1}); % subject number

% customizable options:
runbgncalib = false;
runmidcalib = false;
runendcalib = false;
start_blck  = 1;

if strcmp(argindlg{2},'no')
    start_blck = 1;
end

%% create data folder for subject
foldname = sprintf('./Data/S%02d',subj);
if ~exist(foldname,'dir')
    mkdir(foldname);
end

diary(sprintf('./Data/S%02d/S%02d_Log_%s.txt',subj,subj,datestr(now,'yyyymmdd-HHMM')));

%%
fprintf('*** BEGINNING CALIBRATION ***\n\n')
fprintf('load calibration file...')
[calibfile,calibpath] = uigetfile('*.mat','Select a calibration file');
calibration = importdata(fullfile(calibpath,calibfile));
fprintf('done!\n\n')

if isempty(calibration)
    error('no calibration available');
elseif runbgncalib
    fprintf('Press any key to continue.\n')
    pause
    [calibration1,aborted,errmsg] = run_calibration(subj,calibration,true); % with training!
elseif ~runbgncalib && ~isfield(calibration, 'rslt') 
        fprintf('No calibration results available, run calibration!\n')
        fprintf('Press any key to continue.\n')
        pause
        [calibration1,aborted,errmsg] = run_calibration(subj,calibration,true);
elseif ~runbgncalib && isfield(calibration, 'rslt')
        calibration1 = calibration;   
end
fprintf('\t min / max proportion of Color 1: \t%d / %d\n\n',round(calibration1.rslt.range(1)*100),round(100*calibration1.rslt.range(2)))



fprintf('*** EXPERIMENT 1st Half ***\n\n')
% if strcmp(argindlg{2},'no') || (~strcmp(argindlg{2},'no')&& str2num(argindlg{2})<=8)
% run experiment (1st half or start at specified block)
if ~strcmp(argindlg{2},'no')
    fprintf('load raw experiment structure...')
    [FileName,PathName] = uigetfile('*.mat','Select a raw experiment structure');
    expe_raw = importdata(fullfile(PathName,FileName));
    if isempty(expe_raw)
        error('no raw experiment structure available');
    end
    fprintf('done!\n\n')
else
    expe  = gen_expe(subj);
    shape = shuffle_shape(expe); % create subject shape combinations
    expe_raw       = [];
    expe_raw.expe  = expe;
    expe_raw.shape = shape;
    % save raw expe structure in case of premature termination
    fname = sprintf('DOTCAT_S%02d_%s_raw',subj,datestr(now,'yyyymmdd-HHMM'));
    fname = fullfile(foldname,fname);
    save([fname,'.mat'],'expe_raw');
end
if start_blck <=8
    fprintf('Press any key to continue.\n')
    pause
    [expe,aborted] = run_expe(subj,calibration1,expe_raw,start_blck,8);
    start_blck = 9;
end

fprintf('*** MID-CALIBRATION ***\n\n')
if isempty(calibration)
    error('no calibration available');
elseif runmidcalib
    fprintf('Press any key to continue.\n')
    pause
    [calibration2,aborted,errmsg] = run_calibration(subj,calibration,false);
elseif ~runmidcalib && ~isfield(calibration, 'rslt')
        fprintf('No calibration results available, run calibration!\n')
        fprintf('Press any key to continue.\n')
        pause
        [calibration2,aborted,errmsg] = run_calibration(subj,calibration,false);
elseif ~runmidcalib && isfield(calibration, 'rslt')
        calibration2 = calibration;   
end
fprintf('\t min / max proportion of Color 1: \t%d / %d\n\n',round(calibration2.rslt.range(1)*100),round(100*calibration2.rslt.range(2)))

fprintf('*** EXPERIMENT 2nd Half ***\n\n')
% run experiment (2nd half or start at specified block)
fprintf('Press any key to continue.\n')
pause
[expe,aborted] = run_expe(subj,calibration2,expe_raw,start_blck,12);

fprintf('*** END CALIBRATION ***\n\n')
% run end calibration
if isempty(calibration)
    error('no calibration available');
elseif runendcalib
    fprintf('Press any key to continue.\n')
    pause
    [calibration3,aborted,errmsg] = run_calibration(subj,calibration,false); % no training
elseif ~runendcalib && ~isfield(calibration, 'rslt')
        fprintf('No calibration results available, run calibration!\n')
        fprintf('Press any key to continue.\n')
        pause
        [calibration3,aborted,errmsg] = run_calibration(subj,calibration,false);
elseif ~runendcalib && isfield(calibration, 'rslt')
        calibration3 = calibration;   
end

%
plot_calib(calibration);
savefig(sprintf('./Data/S%02d/calibration1',subj))
plot_calib(calibration2);
savefig(sprintf('./Data/S%02d/calibration2',subj))
plot_calib(calibration3);
savefig(sprintf('./Data/S%02d/calibration3',subj))
%
plot_calib(calibration);
plot_calib(calibration2,true);
plot_calib(calibration3,true);
savefig(sprintf('./Data/S%02d/calibration_both',subj))

fprintf('\nBeginning Calibration:\n')
fprintf('\t min / max proportion of Color 1: \t%d / %d\n\n',round(calibration.rslt.range(1)*100),round(100*calibration.rslt.range(2)))
fprintf('Mid Calibration:\n')
fprintf('\t  min / max proportion of Color 1: \t%d / %d\n\n',round(100*calibration2.rslt.range(1)),round(100*calibration2.rslt.range(2)))
fprintf('End Calibration:\n')
fprintf('\t  min / max proportion of Color 1: \t%d / %d\n\n',round(100*calibration3.rslt.range(1)),round(100*calibration3.rslt.range(2)))



diary('off')