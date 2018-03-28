% Clear the workspace and the screen
close all;
clearvars;
sca

addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/Visual

ppd  = 40;
% set list of color-wise R/G/B values
color_rgb = [ ...
    % http://www.workwithcolor.com/hsl-color-picker-01.htm?
    255,170,128; ... % orange L:75%; Lum:75%
    128,191,255; ... % bleu   L:75%; Lum:75%
    %      255,160,65; ...  % orange AGOBS 255 179 102
    %     65,160,255;...   % blue AGOBS
    136,236,136; ... % green L 73%, Lum 80%, 120°
    253,192,134; ... % orange clair
    136,220,236; ... % blue  L 73%, Lum 80%, 190°
    %255,255,000; ... % yellow
    %255,130,000; ... % orange
    255,000,000; ... % red
    170,000,170; ... % purple
    000,000,255; ... % blue
    035,220,000; ... % green
    255,255,000; ... % yellow
    255,130,000; ... % orange
    ]/255;
color_opa = 2/3; % color opacity
color_shape = [175,175,175]/255;

% set stimulation parameters
lumibg    = 128/255; % background luminance (grey)
fixtn_siz = 0.2*ppd; % fixation point size
hand_siz  = 3*ppd;        % handful size

% set stimulation parameters
shape_siz = 6.0*ppd;        % shape size
shape_off = 7*ppd;          % shape offset
shape_add = round(0.1*ppd); % choice rectangle offset

% set instr parameters
vert_off  = round(1.2*ppd); % vertical offset for instructions screen
instr_shape = 3.0*ppd;
bag_siz = round(6.5*ppd); 
abag_siz = round(.5*bag_siz);
instr_siz = 3.5*ppd;
instr_fac = 1.15;
txtsiz = round(1.0*ppd);

label_agent = 'Piochez un maximum de fois dans le sac';

% set color opacity
color_rgb   = color_rgb*color_opa+lumibg*(1-color_opa);

% set screen parameters
iscr = 0; % screen index
res  = []; % screen resolution
fps  = []; % screen refresh rate

PsychDefaultSetup(2);
video.i = iscr;
video.res = Screen('Resolution',video.i);
video.h = PsychImaging('OpenWindow',video.i,0);
[video.x,video.y] = Screen('WindowSize',video.h);
ifi = Screen('GetFlipInterval', video.h);
Screen('BlendFunction', video.h, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

KbName('UnifyKeyNames');
keywait = KbName('DownArrow');

    shape_tex = zeros(2,8); % shape_tex(1,i)=black (contour) shape_tex(2,i)=greyish (inside)
    instr_tex = zeros(2,8);
    for i = 1:8
        imgc = double(imread(sprintf('./img/shape%dc.png',i)))/255;
        imgc2 = imresize(imgc,instr_siz/size(imgc,1));
        imgc = imresize(imgc,shape_siz/size(imgc,1));
        shape_tex(1,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgc)),imgc),[],[],2);
        instr_tex(1,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgc2)),imgc2),[],[],2);
        imgi = double(imread(sprintf('./img/shape%d.png',i)))/255;
        imgi2 = imresize(imgi,instr_siz/size(imgi,1));
        imgi = imresize(imgi,shape_siz/size(imgi,1));
        shape_tex(2,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgi)),imgi),[],[],2);
        instr_tex(2,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgi2)),imgi2),[],[],2);
    end


% create bag textures 1:Obs | 2:Agent | 3:Agent (in text)
[~,~,imgbagi] = imread('./img/bag6.png'); % filled bag
[~,~,imgbagc] = imread('./img/bagc6.png'); % contour bag
% big bags
imgbagi = double(imgbagi)/255;
imgbagi = imresize(imgbagi,bag_siz/size(imgbagi,1));
imgbagc = double(imgbagc)/255;
imgbagc = imresize(imgbagc,bag_siz/size(imgbagc,1));
bag_tex(1,1) = Screen('MakeTexture',video.h,cat(3,ones(size(imgbagc)),imgbagc),[],[],2);
bag_tex(1,2) = Screen('MakeTexture',video.h,cat(3,ones(size(imgbagi)),imgbagi),[],[],2);

imgbagi = imresize(imgbagi,abag_siz/size(imgbagi,1));
imgbagc = imresize(imgbagc,abag_siz/size(imgbagc,1));
bag_tex(2,1) = Screen('MakeTexture',video.h,cat(3,ones(size(imgbagc)),imgbagc),[],[],2);
bag_tex(2,2) = Screen('MakeTexture',video.h,cat(3,ones(size(imgbagi)),imgbagi),[],[],2);

%create rects
% upper text
Screen('TextSize',video.h,round(txtsiz*instr_fac));
rec_agent = CenterRectOnPoint(Screen('TextBounds',video.h,label_agent),video.x/2-.5*bag_siz,video.y/7);
% bag
bag_rec(1,:) = CenterRectOnPoint(Screen('Rect',bag_tex(1,1)),video.x/2,video.y/2+vert_off-shape_off/2);
bag_rec(2,:) = CenterRectOnPoint(Screen('Rect',bag_tex(1,1)),video.x/2,video.y/2+vert_off+shape_off/2);
bag_rec(3,:) = CenterRectOnPoint(Screen('Rect',bag_tex(2,1)),rec_agent(3)+abag_siz,video.y/7);
% instructions shapes inside of bag
instr_rec(1,:) = CenterRectOnPoint(Screen('Rect',instr_tex(1)),video.x/2,bag_rec(1,2)+.65*bag_siz);
instr_rec(2,:) = CenterRectOnPoint(Screen('Rect',instr_tex(1)),video.x/2,bag_rec(2,2)+.65*bag_siz);

% Over-write the screen in grey
Screen('FillRect', video.h, lumibg);
Screen('Flip',video.h);
WaitKeyPress(keywait);

%Screen('DrawTexture',video.h,bag_tex(1,2),[],bag_rec(1,:),[],[],[],color_rgb(1,:));
Screen('DrawTexture',video.h,bag_tex(1,1),[],bag_rec(1,:),[],[],[],0);
%Screen('DrawTexture',video.h,bag_tex(1,2),[],bag_rec(2,:),[],[],[],color_rgb(2,:));
Screen('DrawTexture',video.h,bag_tex(1,1),[],bag_rec(2,:),[],[],[],0);
% Screen('DrawTexture',video.h,instr_tex(1,2),[],instr_rec(1,:),[],[],[],0);
% Screen('DrawTexture',video.h,instr_tex(2,2),[],instr_rec(1,:),[],[],[],color_shape);
% Screen('DrawTexture',video.h,instr_tex(1,3),[],instr_rec(2,:),[],[],[],0);
% Screen('DrawTexture',video.h,instr_tex(2,3),[],instr_rec(2,:),[],[],[],color_shape);
Screen('DrawTexture',video.h,bag_tex(2,2),[],bag_rec(3,:),[],[],[],color_rgb(1,:));
Screen('DrawTexture',video.h,bag_tex(2,1),[],bag_rec(3,:),[],[],[],0);
Screen('DrawingFinished',video.h);
t = Screen('Flip',video.h);
imgscreen = Screen('GetImage',video.h);
imwrite(imgscreen,'./scrshots/empty_bags.png');

                
WaitKeyPress(keywait);

KbStrokeWait;
sca;
