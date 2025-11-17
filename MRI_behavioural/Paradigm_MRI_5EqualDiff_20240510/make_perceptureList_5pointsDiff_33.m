load('trialStructure.mat')
blockn = 22;

while 1

while 1
%------------------------------------------------------------
% set the distribution of dev for first 3 trials in each block
%------------------------------------------------------------
% basic = Shuffle([zeros(1,blockn);...
%     Shuffle([ones(1,blockn/2),ones(1,blockn/2)*(-1)]);...
%     Shuffle([ones(1,blockn/2)*2,ones(1,blockn/2)*(-2)]);
%     Shuffle([ones(1,blockn/2)*3,ones(1,blockn/2)*(-3)]);
%     Shuffle([ones(1,blockn/2)*4,ones(1,blockn/2)*(-4)])]);

basic = Shuffle([zeros(1,blockn);...%****** changed by xu 20240225
    Shuffle([ones(1,blockn/2),ones(1,blockn/2)*(-1)]);...%****** changed by xu 20240225
    Shuffle([ones(1,blockn/2)*2,ones(1,blockn/2)*(-2)])]);%****** changed by xu 20240225

flag = sum(basic>0);
% fstOneAfterSwtich = abs(basic(1,:));%****** changed by xu 20240323
% if sum(find(flag<1 | flag>3))==0 && sum(fstOneAfterSwtich==1)==sum(fstOneAfterSwtich==2) && sum(fstOneAfterSwtich==2)==sum(fstOneAfterSwtich==3) && sum(fstOneAfterSwtich==3)==sum(fstOneAfterSwtich==4) && sum(fstOneAfterSwtich==4)==sum(fstOneAfterSwtich==0)-2
%     break
% end %****** changed by xu 20240323

fstOneAfterSwtich = (basic(1,:));%****** changed by xu 20240323
if sum(find(flag<0 | flag>2))==0 && sum(fstOneAfterSwtich==2)==sum(fstOneAfterSwtich==1)&& sum(fstOneAfterSwtich==1)== sum(fstOneAfterSwtich==-1)&& ...
    sum(fstOneAfterSwtich==-1)== sum(fstOneAfterSwtich==-2)&& sum(fstOneAfterSwtich==-2)==(sum(fstOneAfterSwtich==0)-2) %****** changed by xu 20240323
    break
end

end

%-----------------------------------------------------------------------------
% set the distribution of dev for last (trialLength(bn)-3) trials in each block
%------------------------------------------------------------------------------
% last = Shuffle([Shuffle([ones(1,5),ones(1,5)*(-1)]);...
%     Shuffle([ones(1,5)*2,ones(1,5)*(-2)]);
%     Shuffle([ones(1,5)*3,ones(1,5)*(-3)]);
%     Shuffle([ones(1,5)*4,ones(1,5)*(-4)])]);

last = Shuffle([Shuffle([ones(1,14),ones(1,14)*(-1)]);...%****** changed by xu 20240225
    Shuffle([ones(1,14)*2,ones(1,14)*(-2)]);
    Shuffle([ones(1,14)*0,ones(1,14)*(0)]);]);%****** changed by xu 20240225

last = last(:);

List   = nan(150,1);
tbasic = 0;
tlast  = 0;
for bn = 1:blockn
    resd = trialLength(bn)-3;%****** changed by xu 20240225
    if resd>0
        List(tbasic+1:tbasic+trialLength(bn)) = [basic(:,bn);Shuffle(last(tlast+1:tlast+resd))];
    else
        List(tbasic+1:tbasic+trialLength(bn)) = [basic(:,bn)];
    end
    tbasic = tbasic+trialLength(bn);
    tlast  = tlast+resd;
end

%---------------------------------------------------------------------------------------
% make rule switch of red to blue: number (3-1, 3-5)=2; blue to red: number (3-1, 3-5)=2
% make rule switch of red to blue and blue to red, without: 1-1, 5-5
%---------------------------------------------------------------------------------------
trialIndex=cumsum(trialLength(1:end-1));
trialIndex_0=trialIndex(find(List(trialIndex)==0)); % before rule switch
trialIndex_2=trialIndex(find(List(trialIndex)==2)); 
trialIndex_neg2=trialIndex(find(List(trialIndex)==-2));% for controlling the trials without simple2simple
% if sum(List(trialIndex_0+1)==-2)==2 && sum(List(trialIndex_0+1)==2)==2  && sum(List(trialIndex_2+1)==2)==0 && sum(List(trialIndex_neg2+1)==-2)==0 ...
%         && sum(mod(ruleList(trialIndex(List(trialIndex_0+1)==-2)),2))==1 && sum(mod(ruleList(trialIndex(List(trialIndex_0+1)==2)),2))==1 % thesame number of trials in red to blue & blue to red
% by xu 20240410
if sum(List(trialIndex_0+1)==-2)>2 && sum(List(trialIndex_0+1)==2)>2 && sum(List(trialIndex_0+1)==2)<=4 && sum(List(trialIndex_0+1)==-2)<=4 ...
        && sum(List(trialIndex_0+1)==-1)>2 && sum(List(trialIndex_0+1)==1)>2 && sum(List(trialIndex_0+1)==1)<=4 && sum(List(trialIndex_0+1)==-1)<=4 ...
        && sum(List(trialIndex_2+1)==2)==0 && sum(List(trialIndex_neg2+1)==-2)==0 % thesame number of trials in red to blue & blue to red
    break
end

end
%---------------------------------------------------------------------------------------
% make actual dev, and randomly set the answer under most difficulty 
%---------------------------------------------------------------------------------------
testList    = List+3;%****** changed by xu 20240225
zeroAngList = Shuffle(find(List==0));
clockwise   = double(List<0);
clockwise(zeroAngList(1:end/2)) = 1;

% save perceptureList testList clockwise
%%
testList_instruct  = testList;
clockwise_instruct = clockwise;
testList_inferred  = testList;
clockwise_inferred = clockwise;
save perceptureList4 testList_instruct clockwise_instruct testList_inferred clockwise_inferred
