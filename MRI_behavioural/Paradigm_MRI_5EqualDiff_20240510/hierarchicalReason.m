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

    
    wid       = 200;% 380;% ****** changed by xu 20240225,addjust according to the distance in MRI or romm ******
    dist      = 100;% 150;% ****** changed by xu 20240225,addjust according to the distance in MRI or romm ****** 
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
    delay1v     = [4 5 6]; %[1 1 1];% %1;%****** changed by xu 20240225
    refer_dur  = .2;
    delay2     = .6;
    test_dur   = .2;
    resp_dur   = 1.2; % changed by wenshan 20211228
    fb_dur     = 1;
    delay3v    = [4 5 6]; %[0.8 1 1.2];%%[2,3,4]; % changed by wenshan 20220105 %****** changed by xu 20240225
    interv     = [4 5 6]; %[1 1.5 2];%%[2,4,6]; % changed by wenshan 20220105 %****** changed by xu 20240225
    
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

% ****************************************************************************************************************

  %% ******************************  
    onset_start(1) = GetSecs;
%     Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
%     Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
%     Screen('Flip', w);
%     WaitSecs(ini_fix/2);

    % New trial
    for i = 1:trialn
        if  i ==51 | i==101 % i ==30 | i==65 | i==96 | i==126  %  ******
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
            Screen('DrawText',w,text_end,HCenter-60,VCenter+60, [255,255,255]);

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
%             Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
%             Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
%             Screen('Flip', w);
%             WaitSecs(ini_fix);
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