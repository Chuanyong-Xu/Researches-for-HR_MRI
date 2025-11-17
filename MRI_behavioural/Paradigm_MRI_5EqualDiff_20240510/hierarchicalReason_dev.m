function hierarchicalReason(sub_name,ID,gender,expIdx,pct,isExp)
% expIdx: 0: instructed; 1: infered
% isExp:  0: practice;1 experiment
% wenshan 20211224
try
    clc;
    Screen('CloseAll');
    path = pwd;
    resultdir  = [pwd,'\result'];
    background = 255;
    %% Pre setting
    % Key
    KbName('UnifyKeyNames')
    EscKey   = KbName('ESCAPE');
    StartKey = KbName('s');
%     Key1     = KbName('2@'); % h
%     Key2     = KbName('1!'); % b
    Key1     = KbName('2@'); % h
    Key2     = KbName('1!'); % b
    %% Setup
    PsychDefaultSetup(1);
    if 1; Screen('Preference', 'SkipSyncTests', 1); end;
    Screen('Preference','TextEncodingLocale','UTF-8');
    HideCursor;sca
    AssertOpenGL;
    scrnNum = max(Screen('Screens'));
    % Screen('Preference', 'ConserveVRAM', 64);
    [w,wRect] = Screen('OpenWindow',scrnNum,0);
    % interval=Screen('GetFlipInterval',w)
    [HCenter,VCenter] = WindowCenter(w);%[wid/2 height/2]
    pixelSize=Screen('PixelSize', w);
    oldTextSize=Screen('TextSize',w,30);

    
    wid       = 200;%380;%****** changed by xu 20240225,addjust according to the distance in MRI or romm ******
    dist      = 100;%150;%****** changed by xu 20240225,addjust according to the distance in MRI or romm ****** 
    fixdis    = 15;
%     X         = ones(100,1)*[-50:49];
%     Y         = [-50:49]'*ones(1,100);
    X         = ones(300,1)*[-150:149];%****** changed by xu 20240402
    Y         = [-150:149]'*ones(1,300);%****** changed by xu 20240402
    Z         = X.^2 + Y.^2;
    
    up_pos    = [HCenter-wid/2 VCenter-dist-wid HCenter+wid/2 VCenter-dist];
    down_pos  = [HCenter-wid/2 VCenter+dist HCenter+wid/2 VCenter+dist+wid];
    centr_pos = [HCenter-wid/2 VCenter-wid/2 HCenter+wid/2 VCenter+wid/2];
    
    cuedir      = cat(2,path,'\Pics');
    instrImage0 = imread([cuedir '\instr_instructed.jpg'],'jpg');
    instrImage1 = imread([cuedir '\instr_inferred.jpg'],'jpg');
    ruleImage0  = imread([cuedir '\blue.jpg'],'jpg');
    ruleImage1  = imread([cuedir '\red.jpg'],'jpg');
    instrIndex0 = Screen('MakeTexture',w,instrImage0); 
    instrIndex1 = Screen('MakeTexture',w,instrImage1);
    ruleIndex0  = Screen('MakeTexture',w,ruleImage0); 
    ruleIndex1  = Screen('MakeTexture',w,ruleImage1);
    
    outcome0     = '\no_reward2.jpg';
    outcome1     = '\coin.jpg';
    targetImage0 = imread([cuedir outcome0],'jpg');
    targetImage1 = imread([cuedir outcome1],'jpg');
    targetIndex0 = Screen('MakeTexture',w,targetImage0);   
    targetIndex1 = Screen('MakeTexture',w,targetImage1);
    %%
    if expIdx == 1
        c = [[255,255,255];[255,255,255];[255,255,255]];% white, white, white
    elseif expIdx == 0
        c = [[0,  0,  255];[255,  0,  0];[255,255,255]]; % blue, red, white
    end
    
     % Trial structure   
    if isExp == 1
        trialn    = 150; %****** changed by xu 20240225
        load('trialStructure.mat')
        rulelist  = ruleList;% red list
        refAngle  = [pct.referAngle{1};pct.referAngle{2};pct.referAngle{3}];
        testAngle = [pct.testedAngle{1};pct.testedAngle{2};pct.testedAngle{3}];
        clockwise = [pct.clockWise{1};pct.clockWise{2};pct.clockWise{3}];
        cList = Shuffle([ones(ceil((trialn-66)*.65),1);zeros(trialn-66-ceil((trialn-66)*.65),1)]);
        t = 0;
        colorList = [];
        for bi = 1:length(trialLength)
            colorList = [colorList;[ones(3,1);cList(t+1:t+trialLength(bi)-3)]];
            t = t+trialLength(bi)-3;
        end
    else
        trialn    = 20; %****** changed by xu 20240225
        rulelist  = double([ones(8,1);ones(4,1)*(-1);ones(8,1)].* (randperm(2,1)-1.5)*2==1);% red list
        t         = randperm(3,1);
        refAngle  = pct.referAngle{t}(end-19:end);
        testAngle = pct.testedAngle{t}(end-19:end);
        clockwise = pct.clockWise{t}(end-19:end);
        colorList = [ones(3,1);Shuffle([ones(2,1);zeros(3,1)]);...
                    ones(3,1);zeros(1,1);...
                    ones(3,1);Shuffle([ones(3,1);zeros(2,1)])];
    end
    answ      = double(rulelist==clockwise);
    up_list   = [zeros(round(trialn/2),1);ones(trialn-round(trialn/2),1)];% red list
    order     = randperm(trialn);
    up_list   = up_list(order);
    
    % color list
    colorlist = nan(trialn,1);
    colorlist(colorList==1) = rulelist(colorList==1)+1;
    colorlist(colorList==0) = 3;
    
    % Duration setting
    ini_fix    = 0; % 10; % initial fixation
    fix_dur    = .5;
    rule_dur   = 1.2; % changed by wenshan 20220105
    delay1v     = [1 1 1];%[4 5 6]; %1;%****** changed by xu 20240225
    refer_dur  = .2;
    delay2     = .6;
    test_dur   = .2;
    resp_dur   = 1.2; % changed by wenshan 20211228
    fb_dur     = 1;
    delay3v    = [0.8 1 1.2];%[4 5 6];%[2,3,4]; % changed by wenshan 20220105 %****** changed by xu 20240225
    interv     = [1 1.5 2];%[4 5 6];%[2,4,6]; % changed by wenshan 20220105 %****** changed by xu 20240225
    trialn_run = round(trialn/3);
    meann      = round(trialn_run/3);
    lown       = round((trialn_run-meann)*.7);
    highn      = trialn_run-meann-lown;
    delay3     = [Shuffle([repmat(delay3v(1),lown,1);repmat(delay3v(2),meann,1);repmat(delay3v(3),highn,1)]);...
                Shuffle([repmat(delay3v(1),lown,1);repmat(delay3v(2),meann,1);repmat(delay3v(3),highn,1)]);...
                Shuffle([repmat(delay3v(1),lown,1);repmat(delay3v(2),meann,1);repmat(delay3v(3),highn,1)])];
    while 1
    iti        = [Shuffle([repmat(interv(1),lown,1);repmat(interv(2),meann,1);repmat(interv(3),highn,1)]);...
                Shuffle([repmat(interv(1),lown,1);repmat(interv(2),meann,1);repmat(interv(3),highn,1)]);...
                Shuffle([repmat(interv(1),lown,1);repmat(interv(2),meann,1);repmat(interv(3),highn,1)])];
            if iti(end)==interv(3);
                break;
            end
    end
    
    delay1   = [Shuffle([repmat(delay1v(1),lown,1);repmat(delay1v(2),meann,1);repmat(delay1v(3),highn,1)]);...
                Shuffle([repmat(delay1v(1),lown,1);repmat(delay1v(2),meann,1);repmat(delay1v(3),highn,1)]);...
                Shuffle([repmat(delay1v(1),lown,1);repmat(delay1v(2),meann,1);repmat(delay1v(3),highn,1)])];%****** changed by xu 20240225

    exp_setting= struct('ini_fix',ini_fix,'fix_dur',fix_dur,'rule_dur',rule_dur,...
        'delay1',delay1,'refer_dur',refer_dur,'delay2',delay2,'test_dur',test_dur,...
        'resp_dur',resp_dur,'fb_dur',fb_dur,'interv',interv,'meann',meann,'lown',lown,...
        'highn',highn,'delay3v',delay3v,'delay3',delay3,'iti',iti,'rulelist',rulelist,...
        'refAngle',refAngle,'testAngle',testAngle,'clockwise',clockwise,'answ',answ,...
        'up_list',up_list,'pct',pct);
    
    onset_start  = nan(3,1);
    onset_fix    = nan(trialn,1);
    onset_rule   = nan(trialn,1);
    onset_delay1 = nan(trialn,1);
    onset_refer  = nan(trialn,1);
    onset_delay2 = nan(trialn,1);
    onset_test   = nan(trialn,1);
    onset_resp   = nan(trialn,1);
    onset_delay3 = nan(trialn,1);
    onset_fb     = nan(trialn,1);
    onset_iti    = nan(trialn,1);
    ruleResp     = nan(trialn,1);
    ruleRT       = nan(trialn,1);
    pcptResp     = nan(trialn,1);
    pcptRT       = nan(trialn,1);

    %%
    % Instruction
    HideCursor;
    tstart = tic;
    ListenChar(0);
    
    % Wait
    eval(['Screen(''DrawTexture'', w, instrIndex' num2str(expIdx) ');']);
    Screen('Flip', w);
    while 1
        [keyIsDown,secs,keyCode]=KbCheck;
        if keyCode(StartKey)||keyCode(EscKey)
            break;
        end
    end
    
    if keyCode(EscKey)
        Screen('CloseAll'); % sca
        cd (path);
        ShowCursor;
        ListenChar(0);
    end
    
%% ****************************** for discrimate the different difficulties ************************
if isExp == 0
        ref_list=[63 64 65 66 67 113 114 115 116 117];%[50:2:60];%********************************
        refAngle_dev  =repmat(ref_list,1,6);
        testAngle_dev( (1-1)*length(ref_list)+1:1*length(ref_list) ) = refAngle_dev( (1-1)*length(ref_list)+1:1*length(ref_list) )+pct.thresh(1);
        testAngle_dev( (2-1)*length(ref_list)+1:2*length(ref_list) ) = refAngle_dev( (2-1)*length(ref_list)+1:2*length(ref_list) )+pct.thresh(2);
        testAngle_dev( (3-1)*length(ref_list)+1:3*length(ref_list) ) = refAngle_dev((3-1)*length(ref_list)+1:3*length(ref_list))+pct.thresh(3);testAngle_dev( (4-1)*length(ref_list)+1:4*length(ref_list) ) = refAngle_dev( (4-1)*length(ref_list)+1:4*length(ref_list) )+pct.thresh(3)
        testAngle_dev( (5-1)*length(ref_list)+1:5*length(ref_list) ) =refAngle_dev((5-1)*length(ref_list)+1:5*length(ref_list))+pct.thresh(4);
        testAngle_dev( (6-1)*length(ref_list)+1:6*length(ref_list) ) = refAngle_dev((6-1)*length(ref_list)+1:6*length(ref_list))+pct.thresh(5);
%         if pct.thresh(3)<0
%             answ_dev(1,:)=[ones(1,length(ref_list)),ones(1,length(ref_list)),Shuffle([ones(1,length(ref_list)),ones(1,length(ref_list))]),zeros(1,length(ref_list)),zeros(1,length(ref_list))];
%             answ_dev(2,:)=[zeros(1,length(ref_list)),zeros(1,length(ref_list)),Shuffle([zeros(1,length(ref_list)),zeros(1,length(ref_list))]),ones(1,length(ref_list)),ones(1,length(ref_list))];
%         elseif pct.thresh(3)>0
%             answ_dev(1,:)=[ones(1,length(ref_list)),ones(1,length(ref_list)),Shuffle([zeros(1,length(ref_list)),zeros(1,length(ref_list))]),zeros(1,length(ref_list)),zeros(1,length(ref_list))];
%             answ_dev(2,:)=[zeros(1,length(ref_list)),zeros(1,length(ref_list)),Shuffle([ones(1,length(ref_list)),ones(1,length(ref_list))]),ones(1,length(ref_list)),ones(1,length(ref_list))];
%         elseif pct.thresh(3)==0
            answ_dev(1,:)=[ones(1,length(ref_list)),ones(1,length(ref_list)),Shuffle(Shuffle([ones(1,length(ref_list)),zeros(1,length(ref_list))])),zeros(1,length(ref_list)),zeros(1,length(ref_list))];
            answ_dev(2,:)=[zeros(1,length(ref_list)),zeros(1,length(ref_list)),Shuffle(Shuffle([zeros(1,length(ref_list)),ones(1,length(ref_list))])),ones(1,length(ref_list)),ones(1,length(ref_list))];
%         end

%         answ_dev(3,:)=[ones(1,6),ones(1,6)*2,ones(1,6)*3,ones(1,6)*4,ones(1,6)*5];
        answ_dev(3,:)=[ones(1,length(ref_list)),ones(1,length(ref_list))*2,ones(1,length(ref_list)*2)*3,ones(1,length(ref_list))*2,ones(1,length(ref_list))*1];
        % random the sequence of the difficulties of different angel 
        angel_answ_all=Shuffle([refAngle_dev;testAngle_dev;answ_dev],1);
        refAngle_dev=angel_answ_all(1,:);
        testAngle_dev=angel_answ_all(2,:);
        answ_dev(1:3,:)=angel_answ_all(3:5,:);

% ****************************************************************************************************************
for dev_rul_i=1:2
        pcptResp_dev=nan(1,length(answ_dev))
        pcptRT_dev=nan(1,length(answ_dev))
    
    text_dev_star = double([['本实验根据您前一天的测试，设定不同难度，']]);
    Screen('DrawText',w,text_dev_star,HCenter-60,VCenter-65, [255,255,255]);
    text_dev_star = double([['现在进行光栅旋转判断练习，请注意区分！']]);
    Screen('DrawText',w,text_dev_star,HCenter-60,VCenter-5, [255,255,255]);
    if dev_rul_i==1
    text_dev_star = double([['假定当前为红色规则：']]);
    Screen('DrawText',w,text_dev_star,HCenter-60,VCenter+65, [255,0,0]);
    text_dev_star = double([['顺时针按“1”键；逆时针按“2”键']]);
    Screen('DrawText',w,text_dev_star,HCenter-60,VCenter+125, [255,0,0]);
    else
    text_dev_star = double([['假定当前为蓝色规则：']]);
    Screen('DrawText',w,text_dev_star,HCenter-60,VCenter+65, [0,0,255]);
    text_dev_star = double([['顺时针按“2”键；逆时针按“1”键']]);
    Screen('DrawText',w,text_dev_star,HCenter-60,VCenter+125, [0,0,255]);
    end
    Screen('Flip', w);
    WaitSecs(5);% WaitSecs(ini_fix/2)

    for dev_i=1:length(testAngle_dev)
        % dev cue
        if answ_dev(3,dev_i)==1
        text_dev_cue = double([['低难度']]);
        Screen('DrawText',w,text_dev_cue,HCenter-4*fixdis, VCenter-fixdis, [255,255,255]);
        elseif answ_dev(3,dev_i)==2
        text_dev_cue = double([['中难度']]);
        Screen('DrawText',w,text_dev_cue,HCenter-4*fixdis, VCenter-fixdis, [255,255,255]);
        elseif answ_dev(3,dev_i)==3
        text_dev_cue = double([['高难度']]);     
        Screen('DrawText',w,text_dev_cue,HCenter-4*fixdis, VCenter-fixdis, [255,255,255]);
        end
        Screen('Flip', w);
        WaitSecs(1); 

        % fix1
        Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        WaitSecs(1);
     % Reference stimuli        
        image = sinGrating([300 300], 2, refAngle_dev(dev_i), 0, 0.5, 0.5);
        image = im2uint8(image);
        image(find(Z >= 150^2)) = 0;
        refImage = Screen('MakeTexture',w,image);
        Screen('DrawTexture', w, refImage,[],centr_pos);
        Screen('Flip', w);
        WaitSecs(refer_dur);
        % fix2
        Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        WaitSecs(delay2);
        % Test stimuli
        image = sinGrating([300 300], 2, testAngle_dev(dev_i), 0, 0.5, 0.5);
        image = im2uint8(image);
        image(find(Z >= 150^2)) = 0;
        testImage = Screen('MakeTexture',w,image);
        Screen('DrawTexture', w, testImage,[],centr_pos);
        Screen('Flip', w);
        WaitSecs(test_dur);
        % Response
        text_qst = double(['?']);
        Screen('DrawText',w,text_qst,HCenter-fixdis,VCenter-20,[255,255,255]);
        Screen('Flip', w);
        onset_resp(dev_i) = GetSecs;

        touch=0;
        while GetSecs-onset_resp(dev_i)<resp_dur && ~(touch & (keyCode(Key1) | keyCode(Key2)))
            [touch, secs, keyCode] = KbCheck;
            if keyCode(Key1) % change here
                pcptResp_dev(dev_i) = 0;
                pcptRT_dev(dev_i)   = secs-onset_resp(dev_i);
            elseif keyCode(Key2)
                pcptResp_dev(dev_i) = 1;
                pcptRT_dev(dev_i)   = secs-onset_resp(dev_i);
            end
        end
        
        while GetSecs-onset_resp(dev_i)<resp_dur
        end
        %
         % Delay3
        Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        WaitSecs(delay3(1));
        
        % Feedback
        if dev_rul_i==1
            if pcptResp_dev(dev_i)==answ_dev(1,dev_i)
                target_Index=eval(['targetIndex1;']);
            else
                target_Index=eval(['targetIndex0;']);
            end
        else
            if pcptResp_dev(dev_i)==answ_dev(2,dev_i)
                target_Index=eval(['targetIndex1;']);
            else
                target_Index=eval(['targetIndex0;']);
            end
        end
        Screen('DrawTexture', w, target_Index,[],centr_pos);
        Screen('Flip', w);
        WaitSecs(fb_dur);
        
        % ITI
        Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        WaitSecs(iti(1));
        
    end

    sess = floor(dev_i/dev_i);% here, max of dev_i ==30
    idx = (sess-1)*dev_i+[1:dev_i];
    if dev_rul_i==1
        acc = length(find(pcptResp_dev(idx)==answ_dev(1,idx)))/length(pcptResp_dev(idx));
        acc_1=length(find(pcptResp_dev(idx)==answ_dev(1,idx)& answ_dev(3,idx)==1))/sum(answ_dev(3,idx)==1);
        acc_2=length(find(pcptResp_dev(idx)==answ_dev(1,idx)& answ_dev(3,idx)==2))/sum(answ_dev(3,idx)==2);
        acc_3=length(find(pcptResp_dev(idx)==answ_dev(1,idx)& answ_dev(3,idx)==3))/sum(answ_dev(3,idx)==3);
%         acc_4=length(find(pcptResp_dev(idx)==answ_dev(1,idx)& answ_dev(3,idx)==4))/length(answ_dev(3,idx)==4);
%         acc_5=length(find(pcptResp_dev(idx)==answ_dev(1,idx)& answ_dev(3,idx)==5))/length(answ_dev(3,idx)==5);
    else
        acc = length(find(pcptResp_dev(idx)==answ_dev(2,idx)))/length(pcptResp_dev(idx));
        acc_1=length(find(pcptResp_dev(idx)==answ_dev(2,idx)& answ_dev(3,idx)==1))/sum(answ_dev(3,idx)==1);
        acc_2=length(find(pcptResp_dev(idx)==answ_dev(2,idx)& answ_dev(3,idx)==2))/sum(answ_dev(3,idx)==2);
        acc_3=length(find(pcptResp_dev(idx)==answ_dev(2,idx)& answ_dev(3,idx)==3))/sum(answ_dev(3,idx)==3);
    end
    miss= length(unique([find(isnan(pcptResp_dev(idx)))]));
    text_end = double(['你刚才漏选了',num2str(miss),'个，正确率是',sprintf('%.2f',acc*100),'% ！']);
    Screen('DrawText',w,text_end,HCenter-60,VCenter-65,[255,255,255]);
    text_end = double(['低难度的正确率是',sprintf('%.2f',acc_1*100),'% ！']);
    Screen('DrawText',w,text_end,HCenter-60,VCenter-5,[255,255,255]);
    text_end = double(['中难度的正确率是',sprintf('%.2f',acc_2*100),'% ！']);
    Screen('DrawText',w,text_end,HCenter-60,VCenter+55,[255,255,255]);
    text_end = double(['高难度的正确率是',sprintf('%.2f',acc_3*100),'% ！']);
    Screen('DrawText',w,text_end,HCenter-60,VCenter+115,[255,255,255]);
    
    Screen('Flip', w);
    WaitSecs(5);% WaitSecs(ini_fix/2);
    
    exp_setting.dev(dev_rul_i).acc(1)= acc; exp_setting.dev(dev_rul_i).acc(2)= acc_1;
    exp_setting.dev(dev_rul_i).acc(3)= acc_2; exp_setting.dev(dev_rul_i).acc(4)= acc_3;
    exp_setting.dev(dev_rul_i).pcptResp_dev=pcptResp_dev;
end
% save data
    exp_setting.dev_angel_answ_all=angel_answ_all;
    

% ****************************************************************************************************************
    text_dev_end = double([['接下来开始层级推理任务练习，注意按键！']]);
    Screen('DrawText',w,text_dev_end,HCenter-60,VCenter-65, [255,255,255]);
    text_dev_end = double([['准备好了按s键开始！']]);
    Screen('DrawText',w,text_dev_end,HCenter-60,VCenter-5, [255,255,255]);
    Screen('Flip', w);

        while 1
            [keyIsDown,secs,keyCode]=KbCheck;
            if keyCode(StartKey)||keyCode(EscKey)
                break;
            end
        end
            
        if keyCode(EscKey)
            Screen('CloseAll'); % sca
            cd (path);
            ShowCursor;
            ListenChar(0);
        end

else
end
% ****************************************************************************************************************

  %% ******************************  
    onset_start(1) = GetSecs;
    Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
    Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
    Screen('Flip', w);
    WaitSecs(ini_fix/2);

    % New trial
    for i = 1:trialn
        if  i ==51 |i==101% ******??????
            save temp_1
            
            sess = floor(i/50);
            idx = (sess-1)*50+[1:50];
            acc = length(find(ruleResp(idx) == exp_setting.rulelist(idx) & pcptResp(idx)==exp_setting.answ(idx)))/length(pcptResp(idx));
            miss= length(unique([find(isnan(ruleResp(idx)));find(isnan(pcptResp(idx)))]));
            text_end = double(['你刚才漏选了',num2str(miss),'个，正确率是',sprintf('%.2f',acc*100),'% ！']);
            Screen('DrawText',w,text_end,HCenter-60,VCenter-65,[255,255,255]);
            
            text_end = double(['请记住你最后选的规则,角度判断,结果！']);%****** by xu 10140513
            Screen('DrawText',w,text_end,HCenter-60,VCenter-5,[255,255,255]);%****** by xu 10140513            
            
            text_end = double([['请休息20秒！准备好了按s键继续！']]);
            Screen('DrawText',w,text_end,HCenter-60,VCenter-5, [255,255,255]);

            Screen('Flip', w);
            WaitSecs(20);
            
            while 1
                [keyIsDown,secs,keyCode]=KbCheck;
                if keyCode(StartKey)||keyCode(EscKey)
                    break;
                end
            end
            
            if keyCode(EscKey)
                Screen('CloseAll'); % sca
                cd (path);
                ShowCursor;
                ListenChar(0);
            end
            
            onset_start(ceil(i/50)) = GetSecs;
            Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
            Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
            Screen('Flip', w);
            WaitSecs(ini_fix);
        end
        % Fixation
        Screen('DrawLine', w, c(colorlist(i),:), HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, c(colorlist(i),:), HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        onset_fix(i) = GetSecs;
        WaitSecs(fix_dur);
        
        % Rule report
        eval(['Screen(''DrawTexture'', w, ruleIndex' num2str(~up_list(i)) ',[],up_pos);']);
        eval(['Screen(''DrawTexture'', w, ruleIndex' num2str(up_list(i))  ',[],down_pos);']);
        Screen('DrawLine', w, c(colorlist(i),:), HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, c(colorlist(i),:), HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        onset_rule(i) = GetSecs;
        
        touch=0;
        [touch, secs, keyCode] = KbCheck; % change here
        while GetSecs-onset_rule(i)<rule_dur && ~(touch & (keyCode(Key1) | keyCode(Key2) | keyCode(EscKey)))
            [touch, secs, keyCode] = KbCheck;
            %              FlushEvents('keyDown');
            if keyCode(Key1) % change here
                ruleResp(i) = double(up_list(i)==0);
                ruleRT(i)   = secs-onset_rule(i);
            elseif keyCode(Key2)
                ruleResp(i) = double(up_list(i)==1);
                ruleRT(i)   = secs-onset_rule(i);
            elseif keyCode(EscKey)
                Screen('CloseAll'); % sca
                cd(path);
                ShowCursor;
                ListenChar(0);
                break;
            end
        end
        
        while GetSecs-onset_rule(i)<rule_dur
        end
        
        % Delay1
        Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        onset_delay1(i) = GetSecs;
        WaitSecs(delay1(i));%****** changed by xu 20240225
        
        % Reference stimuli
        image = sinGrating([300 300], 2, refAngle(i), 0, 0.5, 0.5);
        image = im2uint8(image);
        image(find(Z >= 150^2)) = 0;
        refImage = Screen('MakeTexture',w,image);
        Screen('DrawTexture', w, refImage,[],centr_pos);
        Screen('Flip', w);
        onset_refer(i) = GetSecs;
        WaitSecs(refer_dur);
        
        % Delay2
        Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        onset_delay2(i) = GetSecs;
        WaitSecs(delay2);
        
        % Test stimuli
        image = sinGrating([300 300], 2, testAngle(i), 0, 0.5, 0.5);
        image = im2uint8(image);
        image(find(Z >= 150^2)) = 0;
        testImage = Screen('MakeTexture',w,image);
        Screen('DrawTexture', w, testImage,[],centr_pos);
        Screen('Flip', w);
        onset_test(i) = GetSecs;
        WaitSecs(test_dur);
        
        % Response
        text_qst = double(['?']);
        Screen('DrawText',w,text_qst,HCenter-fixdis,VCenter-20,[255,255,255]);
        Screen('Flip', w);
        onset_resp(i) = GetSecs;
%         clockwise(i)
        %         aaa{1}(i)
        
        touch=0;
        while GetSecs-onset_resp(i)<resp_dur && ~(touch & (keyCode(Key1) | keyCode(Key2)))
            [touch, secs, keyCode] = KbCheck;
            if keyCode(Key1) % change here
                pcptResp(i) = 0;
                pcptRT(i)   = secs-onset_resp(i);
            elseif keyCode(Key2)
                pcptResp(i) = 1;
                pcptRT(i)   = secs-onset_resp(i);
            end
        end
        
        while GetSecs-onset_resp(i)<resp_dur
        end
        
        % Delay3
        Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        onset_delay3(i) = GetSecs;
        WaitSecs(delay3(i));
        
        % Feedback
        if ruleResp(i) == rulelist(i) & pcptResp(i)==answ(i)
            target_Index=eval(['targetIndex1;']);
        else
            target_Index=eval(['targetIndex0;']);
        end
        onset_fb(i) = GetSecs;
        Screen('DrawTexture', w, target_Index,[],centr_pos);
        Screen('Flip', w);
        onset_fb(i) = GetSecs;
        WaitSecs(fb_dur);
        
        % ITI
        Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        onset_iti(i) = GetSecs;
        WaitSecs(iti(i));
    end
    % Save data
    cd(resultdir);
    time_toc = toc(tstart);
    abc      = fix(clock);
    keyword1  = {'instruct','inferred'};
    if isExp
        res_name = ['subj',num2str(ID), '_' keyword1{expIdx+1} '_exp'...
            '_' date '_' num2str(abc(4)) '_' num2str(abc(5)) '_' num2str(abc(6))];
    else
        res_name = ['subj',num2str(ID), '_' keyword1{expIdx+1} '_practice' ...
            '_' date '_' num2str(abc(4)) '_' num2str(abc(5)) '_' num2str(abc(6))];
    end
    eval(['save ' res_name ' time_toc sub_name ID gender trialn exp_setting '...
        'onset_start onset_fix onset_rule onset_delay1 onset_refer onset_delay2 ',...
        'onset_test onset_resp onset_delay3 onset_fb onset_iti ruleResp ruleRT ',...
        'pcptResp pcptRT expIdx;']);
      
    acc = length(find(ruleResp == exp_setting.rulelist & pcptResp==exp_setting.answ))/length(pcptResp);
    miss= length(unique([find(isnan(ruleResp));find(isnan(pcptResp))]));
    text_end = double(['你总共漏选了',num2str(miss),'个，正确率是',sprintf('%.2f',acc*100),'% ！']);
    Screen('DrawText',w,text_end,HCenter-40,VCenter-65,[255,255,255]);
    
    text_end = double(['辛苦啦！任务结束，请联系主试！']);
    Screen('DrawText',w,text_end,HCenter-40,VCenter-5,[255,255,255]);
    Screen('Flip', w);
    while 1
        [keyIsDown,secs,keyCode]=KbCheck;
        if keyCode(EscKey)
            Screen('CloseAll'); % sca
            cd(path);
            ShowCursor;
            ListenChar(0);
            break;
        end
    end
    
    sprintf('Missing trials: %d\n', length(unique([find(isnan(pcptResp));find(isnan(ruleResp))])))

catch
    psychrethrow(psychlasterror);
    Screen('CloseAll');
    ListenChar(0);
    ShowCursor;
    cd (path);
end

end