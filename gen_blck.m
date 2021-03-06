function [blck] = gen_blck(cfg, isprac)
%  GEN_BLCK  Generate block of DOTCAT experiment
%
%  Usage: [blck] = GEN_BLCK(cfg)
%
%  where the configuration structure cfg should contain the following fields:
%   * taskid - task identifier => 1:Observer or 2:Agent
%   * condtn - task condition  => 1:unicolor Dot or 2: bicolor Dot
%   * epimap - episode mapping => 1:color1=target|starting or 2:color2=target|starting
%
%  Output: BLCK structure containing following fields
%   AS IN INPUT CFG:
%   * taskid - task identifier => 1:Observer or 2:Agent
%   * condtn - task condition  => 1:single Dot or 2: bicolor Dot
%   * epimap - episode mapping => 1:color1=target|starting or 2:color2=target|starting
%   NEW FIELDS:
%   * volatility - block volatility (double)
%   * nepi       - number of episodes
%   * ntrl       - number of trial for this block
%   * p_reward   - reward probability => p_reward in [0,1]
%   * reward_seq - probability of reward at each trial => array(1,ntrl)
%   * switch_seq - reward prob switch => 1: switch or 0: no switch
%   * color_seq  - actual rewarding sequence => 1: epimap color or 0: other color
%   * false_seq  - false positives sequence  => 1: false positive or 0: no

if nargin < 1
    error('Missing configuration structure!');
elseif nargin == 1
    isprac = false;
end

% get configuration parameters
taskid = cfg.taskid; % task identifier => 1:observer or 2:agent
condtn = cfg.condtn; % task condition  => 1:single Dot or 2: handful of Dots

% gentype = 2; % switch sequence generation 1: based on volatility | 2: base on nepi

% set block parameters
if  ~isprac
    davg     = 10;
    nepi     = 12;
    p_reward = 0.80; % 1/5 false positive
else % practice
    davg     = 15;
    nepi     = 2;    
    p_reward = 0.90; % 1/10 false positive
end

ntrl   = davg * nepi;  % number of trials per volatility level
dlim   = [5 ntrl];     % min|max number of trials before reversal
rfalse = 1-p_reward;   % proportion of false positives

% generate episodes (Valentin's function)
b = gen_epi(davg,dlim,nepi);

% create block structure
blck             = b;
blck.taskid      = cfg.taskid;
blck.condtn      = cfg.condtn;
blck.epimap      = cfg.epimap;
blck.ntrl        = ntrl;
blck.p_reward    = p_reward;

blck = orderfields(blck,['taskid';'condtn';'epimap';'ntrl'; fieldnames(b);'p_reward']);

% switch sequence
idx_switch = cumsum(b.xs)+1;
switch_seq = zeros(1,ntrl);
switch_seq(idx_switch(1:end-1)) = 1;

% reward probablity sequence
reward_seq = zeros(1,ntrl);
reward_seq(b.ys==1) = p_reward;
reward_seq(b.ys~=1) = 1-p_reward;

%%
blck.reward_seq   = reward_seq;
blck.switch_seq   = switch_seq;

%false positive sequence
false_seq = gen_ffb(blck);

%actual rewarding sequence
color_seq = b.ys;
color_seq(color_seq == 2) = 0;
colorffb_seq = color_seq;
colorffb_seq(false_seq==1) = 1-colorffb_seq(false_seq==1);


blck.color_seq    = color_seq;
blck.colorffb_seq = colorffb_seq;
blck.false_seq    = false_seq;

end