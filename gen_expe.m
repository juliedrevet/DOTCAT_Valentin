function [expe] = gen_expe(subj, taskeq)
%  GEN_EXPE  Generate DOTCAT experiment
%
%  Usage: [expe] = GEN_EXPE(subj,taskeq)
%
%  where subj is the subject number. Experiment parameters are counter-balanced
%  after groups of 8 subjects. The experiment structure expe contains the block
%  substructure blck required to run the experiment and the configuration
%  cfg for this subject
%
%  The optional parameter taskeq controls the degree of equalization between
%  the two tasks (true for full theoretical equalization between the two tasks 
%  (mirror)). The default is true.


addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/

expe  = struct(); %12 => 4 practice blocks, 8 experimental blocks
cfg   = struct();
nprac = 4;

% check input arguments
if nargin < 2
    % apply full equalization between tasks
    taskeq = true;
elseif nargin < 1
    error('Missing subject number!');
end

%% configuration structure
% cfg.taskid: task identifier => 1:observer or 2:agent
% cfg.condtn: task condition  => 1:simple Dot or 2: bicolor dot
% cfg.epimap: episode mapping => 1:color1=target|starting or 2:color2=target|starting


% task identifier => 1:observer or 2:actor
% 1st bit of subject number
if bitget(subj,1)
    cfg.taskid = [1 1 2 2  1 1 2 2  1 1 2 2];
else
    cfg.taskid = [2 2 1 1  2 2 1 1  2 2 1 1];
end
% task condition => 1:simple dot or 2:bicolor dot
% 2nd bit of subject number
if bitget(subj,2)
    cfg.condtn = [1 2 1 2  1 2 1 2  1 2 1 2];
else
    cfg.condtn = [2 1 2 1  2 1 2 1  2 1 2 1];
end
% episode mapping => 1:color1=target|starting or 2:color2=target|starting
% 3rd bit of subject number
if bitget(subj,3)
    cfg.epimap = [2 1 2 1  2 1 1 2  1 2 2 1];
else
    cfg.epimap = [1 2 1 2  1 2 2 1  2 1 1 2];
end

cfg.nprac = nprac;
expe.cfg  = cfg;

nblck = length(cfg.taskid);

%%
b = 1;
fprintf('Generating practice blocks...');
if taskeq
    % full equalization between tasks O and A
    isub = find(cfg.taskid == cfg.taskid(1)); % block indices of 1st task
    jsub = find(cfg.taskid ~= cfg.taskid(1)); % block indices of 2nd task
    for i = 1:length(isub)
        % generate block for 1st task
        iblck = isub(i);
        cfg_blck = [];
        cfg_blck.taskid = cfg.taskid(iblck);
        cfg_blck.condtn = cfg.condtn(iblck);
        cfg_blck.epimap = cfg.epimap(iblck);
        if (iblck <= nprac)
            blck = gen_blck(cfg_blck,true);
            expe.blck(iblck) = blck;
        else
            fprintf('  * generating episodes for block %d/%d...\n',b-nprac,nblck-nprac);
            blck = gen_blck(cfg_blck);
            expe.blck(iblck) = blck;
        end
        b = b+1;
        % generate training blocks independantly (no mirror)
        if iblck <= nprac
            jblck = jsub(i);
            cfg_blck = [];
            cfg_blck.taskid  = cfg.taskid(jblck);
            cfg_blck.condtn  = cfg.condtn(jblck);
            cfg_blck.epimap  = cfg.epimap(jblck);
            blck = gen_blck(cfg_blck,true);
            expe.blck(jblck) = blck;        
        else
            fprintf('  * generating episodes for block %d/%d...\n',b-nprac,nblck-nprac);
            % copy block information for 2nd task in 2nd half
            mirr = jsub(cfg.condtn(jsub)==cfg.condtn(iblck));
            if (iblck<9)
                jblck = mirr(mirr>8);
            else
                jblck = mirr((mirr<9) & (mirr>nprac));
            end
            %expe(jblck).cfg = cfg;
            %expe(jblck).blck = expe(iblck).blck;
            expe.blck(jblck) = expe.blck(iblck);
            %expe(jblck).blck.taskid = 3-expe(iblck).blck.taskid;
            expe.blck(jblck).taskid = 3-expe.blck(iblck).taskid;
        end
        if b == nprac
            fprintf('done.\nGenerating testing blocks...\n');
        elseif b == length(expe)
            fprintf('...done.\n\n')
        end
        b = b+1;
    end
else
    % no equalization between tasks
    for iblck = 1:nblck
        %expe(iblck).cfg = cfg;
        % generate block
        cfg_blck = [];
        cfg_blck.taskid = cfg.taskid(iblck);
        cfg_blck.condtn = cfg.condtn(iblck);
        cfg_blck.epimap = cfg.epimap(iblck);
        if iblck <= nprac
            %fprintf('  * generating episodes for practice block %d/%d...\n',b,nprac);
            blck = gen_blck(cfg_blck,true);
        else
            fprintf('  * generating episodes for block %d/%d...\n',b-nprac,nblck-nprac);
            blck = gen_blck(cfg_blck);
        end
        blck.taskid = cfg_blck.taskid;
        %expe(iblck).blck = blck;
        expe.blck(iblck) = blck;
        if b == nprac
            fprintf('done.\nGenerating testing blocks...\n');
        elseif b == length(expe)
            fprintf('...done.\n\n')
        end
        b = b+1;
    end
end

%plot_expe(expe)
end