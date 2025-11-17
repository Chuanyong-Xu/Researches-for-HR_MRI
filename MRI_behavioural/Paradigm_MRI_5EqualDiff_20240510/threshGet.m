function [thresh] = threshGet(sub_name,ID,gender, isExp)
% isExp: 0: practice;1 experiment
% wenshan 20211224
try
    % Clear the workspace
    clc;
    Screen('CloseAll');
    
    path = pwd;
    resultdir  = [pwd,'\result'];
    background = 255;
    c = [255,255,255];
    
    if isExp == 1
        trialn     = 300;%360;%328;%****** changed by xu 20240402
        % before setting
        orRange    = 18;%16.5;%5;%21; %****** changed by xu 20240402
        numSteps   = 10;%12;%6;%41; % change %****** changed by xu 20240402
    else
        trialn     = 20;
        % before setting
        orRange    = 30; %
        numSteps   = 5; % change
    end

%----------------------------------------------------------------------
%                       Keyboard information
%----------------------------------------------------------------------

    KbName('UnifyKeyNames')
    EscKey   = KbName('ESCAPE');
    StartKey = KbName('S');
    Key1     = KbName('h'); %change here%*********
    Key2     = KbName('b');%*******
    
%----------------------------------------------------------------------
%                       Timing Information
%----------------------------------------------------------------------
    % Duration setting
    
    ini_fix    = 1; % initial fixation
    fix_dur    = .5;
    refer_dur  = .2;
    delay2     = .6;
    test_dur   = .2;
    resp_dur   = 1.2; % changed by wenshan 20211228
    delay3v     = [0.8 1 1.2];%[4 5 6]; %.5;%****** changed by xu 20240225
    fb_dur     = 1;
    
   % iti
    interv     = [1 1.5 2];%[4 5 6];%[.6,.8,1];%****** changed by xu 20240225
    meann      = round(trialn/3);
    lown       = round((trialn-meann)*.7);
    highn      = trialn-meann-lown;
    timeList   = [repmat(interv(1),lown,1);repmat(interv(2),meann,1);repmat(interv(3),highn,1)];
    iti        = timeList(randperm(trialn));
    %delay3 %****** changed by xu 20240225
    timeList   = [repmat(delay3v(1),lown,1);repmat(delay3v(2),meann,1);repmat(delay3v(3),highn,1)];%****** changed by xu 20240225
    delay3    = timeList(randperm(trialn));%****** changed by xu 20240225
    
    %--------------------
    % Gabor information
    %--------------------
    
    % reference
    ref_angle  = 65;
    addAng     = [-2 -1 0 1 2 48 49 50 51 52]';%[-5:1:5]';%[-1 0 1]';%[-5:1:5]';%****** changed by xu 20240402
    addAngList = Shuffle([repmat(addAng,floor(trialn/length(addAng)),1);...
        zeros(trialn-length(addAng)*floor(trialn/length(addAng)),1)]);
    refAngle   = ref_angle + addAngList;

    % testAng 
%     testAng = linspace(-orRange, orRange, numSteps)';
%     testAngList = round(rand(trialn,1)*orRange*2)-orRange;
    testAng     = linspace(-orRange, orRange, numSteps)';

    testAngList = Shuffle([repmat(testAng,floor(trialn/length(testAng)),1);...
        zeros(trialn-length(testAng)*floor(trialn/length(testAng)),1)]);
%     zeroAngList = Shuffle(find(testAngList==0));%****** changed by xu 20240402
%     histogram(testAngList)
    
%     testAng = linspace(-orRange, orRange, numSteps)';
%     testAngList= Shuffle([repmat(testAng,floor(trialn/length(testAng)),1);zeros(trialn-length(testAng)*floor(trialn/length(testAng)),1)]);

    testAngle  = refAngle + testAngList;
    clockwise  = double(testAngList<0);
% % % % % %     clockwise(zeroAngList(1:end/2)) = 1;%****** changed by xu 20240402

    respVector  = zeros(1, size(testAng,1));
    countVector = zeros(1, size(testAng,1));  
    
    exp_setting= struct('fix_dur',fix_dur,'refer_dur',refer_dur,'delay2',delay2,'test_dur',test_dur,...
        'resp_dur',resp_dur,'fb_dur',fb_dur,'interv',interv,'meann',meann,'lown',lown,...
        'highn',highn,'timeList',timeList,'iti',iti,...
        'refAngle',refAngle,'testAngle',testAngle,'clockwise',clockwise);
    %% Setup
    PsychDefaultSetup(1);
    if 1; Screen('Preference', 'SkipSyncTests', 1); end;
    Screen('Preference','TextEncodingLocale','UTF-8');
    HideCursor;sca
    AssertOpenGL;
    scrnNum = max(Screen('Screens'));
%     Screen('Preference', 'ConserveVRAM', 64);
    [w,wRect] = Screen('OpenWindow',scrnNum,0);
    % interval=Screen('GetFlipInterval',w)
    [HCenter,VCenter] = WindowCenter(w);
    oldTextSize=Screen('TextSize',w,30);
    %% Pre setting
    wid       = 200;
    dist      = 100;
    fixdis    = 15;
%     X         = ones(100,1)*[-50:49];
%     Y         = [-50:49]'*ones(1,100);
    X         = ones(300,1)*[-150:149];%****** changed by xu 20240402
    Y         = [-150:149]'*ones(1,300);%****** changed by xu 20240402
    Z         = X.^2 + Y.^2;
    
    centr_pos= [HCenter-wid/2 VCenter-wid/2 HCenter+wid/2 VCenter+wid/2];
    
    cuedir      = cat(2,path,'\Pics');
    instrImage0 = imread([cuedir '\instr_threshold.jpg'],'jpg');
    instrIndex0 = Screen('MakeTexture',w,instrImage0); 
    
    outcome0 = '\no_reward2.jpg';
    outcome1 = '\coin.jpg';
    targetImage0 = imread([cuedir outcome0],'jpg');
    targetImage1 = imread([cuedir outcome1],'jpg');
    targetIndex0 = Screen('MakeTexture',w,targetImage0);   
    targetIndex1 = Screen('MakeTexture',w,targetImage1);
    
    onset_fix    = nan(trialn,1);
    onset_refer  = nan(trialn,1);
    onset_delay2 = nan(trialn,1);
    onset_test   = nan(trialn,1);
    onset_resp   = nan(trialn,1);
    onset_delay3 = nan(trialn,1);
    onset_fb     = nan(trialn,1);
    onset_iti    = nan(trialn,1);
    pcptResp     = nan(trialn,1);
    pcptRT       = nan(trialn,1);

    %%
    tstart = tic;
    ListenChar(1);
    
    % Wait
    Screen('DrawTexture', w, instrIndex0);
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
    onset_start = GetSecs;
    WaitSecs(ini_fix);

%----------------------------------------------------------------------
%                       Experimental loop
%----------------------------------------------------------------------

    HideCursor;
    for i = 1:trialn
%         if i == 100 | i == 200 %****** changed by xu 20240225
%          if i == 82 | i == 164 | i == 246 %****** changed by xu 20240402
        if i == 120 | i == 240
            text_end = double(['辛苦啦，休息一会，准备好了按s键继续实验！']);
            Screen('DrawText',w,text_end,HCenter-40,VCenter-5,[255,255,255]);%****** changed by xu 20240225
            Screen('Flip', w);
            while 1
                [keyIsDown,secs,keyCode]=KbCheck;
                if keyCode(StartKey)
                    break;
                end
            end
        end
        % Fixation
        Screen('DrawLine', w, c, HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, c, HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        onset_fix(i) = GetSecs;
        WaitSecs(fix_dur);
      
        % Reference stimuli
%         image = sinGrating([100 100], 9, refAngle(i), 0, 0.5, 0.5);%****** changed by xu 20240402
        image = sinGrating([300 300], 2, refAngle(i), 0, 0.5, 0.5); %****** changed by xu 20240402
        image = im2uint8(image);
%         image(find(Z >= 50^2)) = 0; %****** changed by xu 20240402
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
        image = sinGrating([300 300], 2, testAngle(i), 0, 0.5, 0.5); %****** changed by xu 20240402
        image = im2uint8(image);
        image(find(Z >= 150^2)) = 0; %****** changed by xu 20240402
        testImage = Screen('MakeTexture',w,image);
        Screen('DrawTexture', w, testImage,[],centr_pos);
        Screen('Flip', w);
        onset_test(i) = GetSecs;
        WaitSecs(test_dur);
        
        % Response
        text_qst = double('?');
        Screen('DrawText',w,text_qst,HCenter-fixdis,VCenter-20,[255,255,255]);
        Screen('Flip', w);
        onset_resp(i) = GetSecs;
        clockwise(i)
        
        touch=0;
%         [touch, secs, keyCode] = KbCheck; % change here %****** changed by xu 20240225
        while GetSecs-onset_resp(i)<resp_dur && ~(touch & (keyCode(Key1) | keyCode(Key2))) 1 
            [touch, secs, keyCode] = KbCheck;
            if keyCode(Key1) % change here
                pcptResp(i) = 0;
                pcptRT(i)   = secs-onset_resp(i);
            elseif keyCode(Key2)
                pcptResp(i) = 1;
                pcptRT(i)   = secs-onset_resp(i);
            end
        end 
        
        while GetSecs-onset_resp(i)<resp_dur %****** changed by xu 20240225
        end %****** changed by xu 20240225
%         
        % Delay3
        Screen('DrawLine', w, [255,255,255], HCenter-fixdis,VCenter, HCenter+fixdis, VCenter, 5);
        Screen('DrawLine', w, [255,255,255], HCenter, VCenter-fixdis, HCenter, VCenter+fixdis, 5);
        Screen('Flip', w);
        onset_delay3(i) = GetSecs;
        WaitSecs(delay3(i));%****** changed by xu 20240225
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
        % Feedback
        if pcptResp(i) == clockwise(i)
            target_Index = eval(['targetIndex1;']);
        else
            target_Index = eval(['targetIndex0;']);
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
        
        % Record
          respVector(testAng == testAngList(i)) = respVector(testAng == testAngList(i)) + (pcptResp(i)==0);

        % Add one to the counter for that stimulus
        if(~isnan(pcptResp(i)))
            countVector(testAng == testAngList(i)) = countVector(testAng == testAngList(i)) + double(~isnan(pcptResp(i))); %changed by wenshan,20211227
        end
        
    end
      
    
    % --------------------------------------
    % save data
    % --------------------------------------
    cd(resultdir);
    time_toc = toc(tstart);
    abc      = fix(clock);
    keyword  = {'practice','exp'};
    res_name = ['subj',num2str(ID), '_threshGet_',keyword{isExp+1} ...
        '_' date '_' num2str(abc(4)) '_' num2str(abc(5)) '_' num2str(abc(6))];
    eval(['save ' res_name ' time_toc sub_name ID gender trialn exp_setting '...
        'onset_start onset_fix  onset_refer onset_delay2 ',...
        'onset_test onset_resp onset_delay3 onset_fb onset_iti ',...
        'pcptResp pcptRT testAng respVector countVector;']);

    
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

        
    % --------------------------------------
    % get threshold
    % --------------------------------------
%     b = glmfit(testAng, [respVector countVector], 'binomial', 'link', 'probit');
%     yfit = glmval(b, testAng, 'probit', 'binomialsize', countVector);
%     figure;
%     plot(testAng, respVector ./ countVector, 'o', testAng, yfit./countVector, '-')
%     axis([min(testAng(:, 1)) max(testAng(:, 1)) 0 1]);
%     xlabel('Angle of Orientation (Degrees)');
%     ylabel('Performance');
%     title('Psychometric function');
    
%   thrs = [.01, 1/16, 1/8, 1/4, 1/2, 3/4, 7/8, 15/16, .99];
% %     thrs = [1/6, 2/6, 3/6, 4/6, 5/6]; %****** changed by xu 20240225
    thrs= [0.1, 0.3, 0.5, 0.7, 0.9]; %****** changed by xu 20240305
    
    data = [testAng'; respVector; countVector]';
    [co, curve, thresh] = ...
        fitPsyche.fitPsycheCurveLogit(data(:, 1), data(:, 2) ./ data(:, 3),'targets',thrs);
    figure;
%     yfit = glmval(co, data(:, 1), 'logit', 'binomialsize', data(:, 3));
    plot(thresh, thrs,'ro',data(:, 1), data(:, 2) ./ data(:, 3), 'ko')
%     plot(data(:, 1), data(:, 2) ./ data(:, 3), 'ro-', 'MarkerFaceColor', 'r');
    hold on
    plot(curve(:,1),curve(:,2),'MarkerFaceColor', 'r');
    plot([0 0],[0 1])
    hold on;
    plot([min(data(:,1)) max(data(:,1))],[0.5  0.5])
%     box off
%     axis([min(data(:, 1)) max(data(:, 1)) 0 1]);
    axis([min(min(data(:, 1)),thresh(1)) max(max(data(:, 1)),thresh(end)) 0 1]);
    xlabel('Angle of Orientation (Degrees)','fontsize',18);
    ylabel('Performance','fontsize',18);
    title('Psychometric function by GLM','fontsize',20);
 %****************************************************************************% threshold calculated by WH2001
     % UpperLimits:
    UL = [0.01, 0.01,18,18]; % Limit upper bound of g and l to 5%
    %     UL = [0.01, 0.01,16.5,16.5]; % Limit upper bound of g and l to 5%
    % StartPoints:
    SP = [0, 0, 0.5, 0.5];
    % LowerLimits:
    LL = [0.01, 0.01, -18,-18];
%     LL = [0.01, 0.01, -16.5,-16.5];

    ffit1=fitPsyche(data(:, 1), data(:, 2) ./ data(:, 3),'WH',[UL;SP;LL]);
%     ffit1=fitPsyche(data(:, 1), data(:, 2) ./ data(:, 3),'WH');
        display(ffit1.model)
        figure
        plotPsyche(ffit1)
            hold on;
            plot(data(:,1),(data(:, 2) ./ data(:, 3)),'o')
            hold on;
            plot([0 0],[0 1])
            hold on;
            plot([min(data(:,1)) max(data(:,1))],[0.5  0.5])
        legend({'y1', 'y1 fit'}, 'Location', 'NorthWest')
        title('WH 2001 fit')
 %****************************************************************************% threshold calculated by WH2001   
    sprintf('Missing trials: %d\n', sum(isnan(pcptResp)))
    catch
    psychrethrow(psychlasterror);
    Screen('CloseAll');
    ListenChar(0);
    ShowCursor;
end
end
