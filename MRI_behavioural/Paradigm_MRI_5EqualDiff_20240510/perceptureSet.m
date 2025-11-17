function [pct] = perceptureSet(threshold, expIdx, ID)
% threshold: subjective threshold
% expIdx: experimental type. 0: instructed; 1: infered.
% wenshan, 20211222

ref_angle  = 65;%55;
runn       = 3;
% addAng     = [-1 0 1]';%[-5:1:5]';%****** changed by xu 20240402
addAng     = [-2 -1 0 1 2 48 49 50 51 52]';%[-5:1:5]';%****** changed by xu 20240402
trialn     = 50;

% file = ['./pcptList/perceptureList',num2str(mod(ID,5))];% 10 changged to 5 by wenshan 20220224
file = ['./pcptList/perceptureList',num2str(mod(ID,1))];% 5 changged to 1 by xu
load(file)

if expIdx == 0
    testList  = testList_instruct;
    clockwise = clockwise_instruct;
elseif expIdx == 1
    testList  = testList_inferred;
    clockwise = clockwise_inferred;
end

times1     = floor(trialn/length(addAng));
addAngList = [repmat(addAng,times1,1);zeros(trialn-length(addAng)*times1,1)];
testAngList= threshold(testList);

for ri = 1:runn
    addAngList     = addAngList(randperm(trialn));
    pct.referAngle{ri} = ref_angle + addAngList;
    pct.testedAngle{ri}= pct.referAngle{ri} + [testAngList((ri-1)*trialn+1:ri*trialn)]';
    pct.clockWise{ri}  = clockwise((ri-1)*trialn+1:ri*trialn);
end

% save perceptureList referAngle testedAngle clockWise aaa
end


