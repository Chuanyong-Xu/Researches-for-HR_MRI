clear; clc
restoredefaultpath
addpath('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/MRI_behavioural-29subs/Hiearchica_Reasoning_model/depend/Functions');
load('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/DeepRL_RNN/data/update_data_switch_weight_diff0.mat');
%collect data
for subs =1:length(infer_data)
    pr_of_swi{subs,1}  = double(infer_data{subs}.tDev)';
    pr_of_swi{subs,2}  = double(infer_data{subs}.Nback)';
    pr_of_swi{subs,3}  = double(infer_data{subs}.pr_of_switch)';
    pr_of_swi{subs,4}  = double(infer_data{subs}.rnn_pre);
    pr_of_swi{subs,5}  = double(infer_data{subs}.hidden_state);
    for con=1:3
        Pr_Sw_o_3_baye{con}(subs,1) = mean(pr_of_swi{subs,3}(pr_of_swi{subs,1}==1 &pr_of_swi{subs,2}==con-1));
        Pr_Sw_o_3_baye{con}(subs,2) = mean(pr_of_swi{subs,3}(pr_of_swi{subs,1}==2 &pr_of_swi{subs,2}==con-1));
        Pr_Sw_o_3_baye{con}(subs,3) = mean(pr_of_swi{subs,3}(pr_of_swi{subs,1}==3 &pr_of_swi{subs,2}==con-1));

        Pr_Sw_inp_3_rnn{con}(subs,1) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==1 &pr_of_swi{subs,2}==con-1));
        Pr_Sw_inp_3_rnn{con}(subs,2) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==2 &pr_of_swi{subs,2}==con-1));
        Pr_Sw_inp_3_rnn{con}(subs,3) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==3 &pr_of_swi{subs,2}==con-1));
    end

    conds_9(subs,1,:)= mean(pr_of_swi{subs,5}(pr_of_swi{subs,1}==1 &pr_of_swi{subs,2}==0,:));
    conds_9(subs,2,:)= mean(pr_of_swi{subs,5}(pr_of_swi{subs,1}==2 &pr_of_swi{subs,2}==0,:));
    conds_9(subs,3,:)= mean(pr_of_swi{subs,5}(pr_of_swi{subs,1}==3 &pr_of_swi{subs,2}==0,:));
    conds_9(subs,4,:)= mean(pr_of_swi{subs,5}(pr_of_swi{subs,1}==1 &pr_of_swi{subs,2}==1,:));
    conds_9(subs,5,:)= mean(pr_of_swi{subs,5}(pr_of_swi{subs,1}==2 &pr_of_swi{subs,2}==1,:));
    conds_9(subs,6,:)= mean(pr_of_swi{subs,5}(pr_of_swi{subs,1}==3 &pr_of_swi{subs,2}==1,:));
    conds_9(subs,7,:)= mean(pr_of_swi{subs,5}(pr_of_swi{subs,1}==1 &pr_of_swi{subs,2}>=2,:));
    conds_9(subs,8,:)= mean(pr_of_swi{subs,5}(pr_of_swi{subs,1}==2 &pr_of_swi{subs,2}>=2,:));
    conds_9(subs,9,:)= mean(pr_of_swi{subs,5}(pr_of_swi{subs,1}==3 &pr_of_swi{subs,2}>=2,:));    
end

save RNN_midLayer_conds_9.mat conds_9

%% plotting
    set(0, 'DefaultAxesFontName', 'Helvetica');   % Fonttype for axis
    set(0, 'DefaultTextFontName', 'Helvetica');   % Fonttype for text

addpath(genpath('/Users/vincentxu/Desktop/MRI_behavioural/Hiearchica_Reasoning_model/depend/Functions'));

load('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/model_res_switch_fmri_m_4_29sub.mat')
%
for iT = 1:3
    res_o_switch_3 =  [res_o{iT}.res_switch(:,3),res_o{iT}.res_switch(:,2)+res_o{iT}.res_switch(:,4),res_o{iT}.res_switch(:,1)+res_o{iT}.res_switch(:,5)];

    res_o_count_3 = [res_o{iT}.res_count(:,3),res_o{iT}.res_count(:,2)+res_o{iT}.res_count(:,4),res_o{iT}.res_count(:,1)+res_o{iT}.res_count(:,5)];

    Pr_Sw_o_3{iT} = res_o_switch_3 ./ res_o_count_3;%******
end
% figure_1: 3 point
id=1:29;
for ki = 1:3
    Pr_Sw_o_3{ki}       = Pr_Sw_o_3{ki}(1:29,:);
    Pr_Sw_o_3_baye{ki}  = Pr_Sw_o_3_baye{ki}(1:29,:);
    Pr_Sw_inp_3_rnn{ki} = Pr_Sw_inp_3_rnn{ki}(1:29,:);
end
    Pr_Sw_o_3_new{1} = Pr_Sw_o_3{3};
    Pr_Sw_o_3_new{2} = Pr_Sw_o_3{1};
    Pr_Sw_o_3_new{3} = Pr_Sw_o_3{2};
    Pr_Sw_o_3        = Pr_Sw_o_3_new;

c = {[0.4660 0.6740 0.1880], [0.2784    0.4471    0.451], [0.7412    0.3569    0.0235]};%{[0.9290 0.6940 0.1250],[0.6350 0.0780 0.1840]};
legendt = {'1B-Er: Data','2B-Er: Data','1B-Er: Model','2B-Er: Model'};

figure
set(gcf,'unit','pixels','position',[200 50 450 400]);
subplot(1,1,1)
    set(gca,'FontSize',14,'tickdir','out')
for ki = 1:3
    ste = nanstd(Pr_Sw_o_3{ki})./sqrt((length(id)-sum(isnan(Pr_Sw_o_3{ki}))));
    be(ki)=errorbar(1:3, nanmean(Pr_Sw_o_3{ki}),ste, 'ko','LineWidth',2,'MarkerFaceColor',c{ki},'MarkerSize',16); hold on
end
rn(1)=fillsteplottusc(Pr_Sw_inp_3_rnn{2},4); hold on
rn(2)=fillsteplottust(Pr_Sw_inp_3_rnn{3},4); hold on
rn(3)=fillsteplotgreen(Pr_Sw_inp_3_rnn{1},4); hold on

set(gca,'XTiCk',1:1:3);
set(gca,'XTickLabel',{'Hard (D1)' 'Middle (D2)' 'Easy (D3)'});
axis([0.5 3.5 0 1]);yticks([0:.25:1])
ylabel('Pr(isSwitch)','fontsize',14); xlabel('Difficulties','fontsize',14)
set(gca,'FontSize',14,'LineWidth',1)
set(gca,'tickdir','out')
set(gca,'box','off')
legend([be(3),be(2),be(1),rn(2),rn(1),rn(3)],{'E2','E1','E0','E2 RNN','E1 RNN','E0 RNN'},...
    'Location','northeastoutside','Orientation','vertical');  legend boxoff;
title('Participants & RNN','fontsize',18,'FontWeight','Normal')

print('Pr_switch_bayesian_rnn0', '-dsvg','-painters')


%% figure 2
%MDS
addpath(genpath('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/matlab-special-heatmap'));
colorlist = flip(slanCM(188));
%plot RDM----------------------------------------------------------------
data_RDM = reshape(mean(conds_9,1),size(conds_9,2),1*size(conds_9,3));
DisMatrix = squareform(pdist(data_RDM, 'mahal',covdiag(data_RDM)));
% DisMatrix = squareform((pdist(data_RDM, 'correlation')));

addpath(genpath('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN/matlab-special-heatmap'));
figure('visible','on','position',[350,200,400,300])
    imagesc(DisMatrix);
    hold on
%     caxis([0,3]); % for newest matlab version: clim([-1,1]);
    colormap(flipud(slanCM(21))); colorbar %97
    axis([.5 9.5 .5 9.5]);
    xticks(1:1:9); xticklabels({'0er-hard','0er-mid','0er-easy','1er-hard','1er-mid','1er-easy','2er-hard','2er-mid','2er-easy'}); 
    yticks(1:1:9); yticklabels({'0er-hard','0er-mid','0er-easy','1er-hard','1er-mid','1er-easy','2er-hard','2er-mid','2er-easy'});
    ax=gca; ax.FontSize=10; ax.XTickLabelRotation=90; box off
    title('RNN RDM')
    hold off
print(['figures/' 'RNN_correlation_RDM'], '-dsvg','-painters');
%------------------------------------------------------------------------


%% 2D
% data_RDM = reshape((beta_cond_subs),size(beta_cond_subs,2),size(beta_cond_subs,1)*size(beta_cond_subs,3));
data_RDM = reshape(mean(conds_9,1),size(conds_9,2),1*size(conds_9,3));
DisMatrix = squareform(pdist(data_RDM, 'mahal',covdiag(data_RDM)));

% data_RDM=smoothRDM(data_RDM,9) %*************************
% DisMatrix = squareform((pdist(data_RDM, 'euclidean')));

figure; set(gcf,'unit','pixels','position',[500 50 1200 200])
% 3-D
[Y, stress,dispa] = mdscale(DisMatrix,3);
Y=zscore(Y);

% [Y, stress,dispa] = mdscale(DisMatrix,3,'start','random');
% [Y, stress] = cmdscale(DisMatrix);
% [Y, y] = pcm_classicalMDS(DisMatrix);
subplot(1,4,1)
for p_i = 1:size(Y,1)
    s=scatter3(Y(p_i,1),Y(p_i,2),Y(p_i,3),'filled');
    s.SizeData = 300; s.MarkerFaceColor=colorlist(p_i*15,:)
    hold on
    text(Y(p_i,1),Y(p_i,2),Y(p_i,3),sprintf('%.0f %.0f %.0f',round(p_i)))
end
% axis equal
% axis([-1 1 -1 1 -1 1])
fill3([Y(1,1),Y(2,1),Y(3,1)],...
      [Y(1,2),Y(2,2),Y(3,2)],...
      [Y(1,3),Y(2,3),Y(3,3)],colorlist(90,:),'EdgeColor','none')
fill3([Y(4,1),Y(5,1),Y(6,1)],...
      [Y(4,2),Y(5,2),Y(6,2)],...
      [Y(4,3),Y(5,3),Y(6,3)],colorlist(130,:),'EdgeColor','none')
fill3([Y(7,1),Y(8,1),Y(9,1)],...
      [Y(7,2),Y(8,2),Y(9,2)],...
      [Y(7,3),Y(8,3),Y(9,3)],colorlist(160,:),'EdgeColor','none') 
grid on
% legend('1:rew-Hard','2:rew-Mid','3:rew-Easy','4:1E-Hard','5:1E-Mid','6:1E-Easy','7:2E-Hard','8:2E-Mid','9:2E-Easy')

xlabel('MDS1'); ylabel('MDS2'); zlabel('MDS3');
% view(40,30)
hold off

% 3-D
subplot(1,4,2)
for p_i = 1:size(Y,1)
    s=scatter(Y(p_i,1),Y(p_i,2),'filled');
    s.SizeData = 300; s.MarkerFaceColor=colorlist(p_i*15,:)
    hold on
    text(Y(p_i,1),Y(p_i,2),sprintf('%.0f %.0f',round(p_i)))
end
% axis equal
grid on
% legend('Hard','Mid','Easy','0Er','1E','2E','Diff','Ners')
xlabel('MDS1'); ylabel('MDS2');
hold off

% 3-D
subplot(1,4,3)
for p_i = 1:size(Y,1)
    s=scatter(Y(p_i,1),Y(p_i,3),'filled');
    s.SizeData = 300; s.MarkerFaceColor=colorlist(p_i*15,:)
    hold on
    text(Y(p_i,1),Y(p_i,3),sprintf('%.0f %.0f',round(p_i)))
end
    plot(Y([2 1 3],1),Y([2 1 3],3),'-','LineWidth',2,'Color',colorlist(30,:));
    plot(Y([5 4 6],1),Y([5 4 6],3),'-','LineWidth',2,'Color',colorlist(90,:));
    plot(Y([8 7 9],1),Y([8 7 9],3),'-','LineWidth',2','Color',colorlist(130,:));
    %for 0 1 2 errs
    plot(Y([1 4 7],1),Y([1 4 7],3),'-','LineWidth',3','Color',colorlist(130,:));
grid on
xlabel('MDS1'); ylabel('MDS3');
hold off

% 3-D
subplot(1,4,4)
for p_i = 1:size(Y,1)
    s=scatter(Y(p_i,2),Y(p_i,3),'filled');
    s.SizeData = 300; s.MarkerFaceColor=colorlist(p_i*15,:)
    hold on
    text(Y(p_i,2),Y(p_i,3),sprintf('%.0f %.0f',round(p_i)))
end
grid on
xlabel('MDS2'); ylabel('MDS3');
hold off

print(['GLM2_Dev_Ners_Geome_MDS_RNN0'], '-dsvg','-painters');





%% figure 2 for PCA
%MDS
addpath(genpath('/Users/vincentxu/Desktop/RNN/matlab-special-heatmap'));
colorlist = flip(slanCM(188));
%2D
% data_RDM = reshape((beta_cond_subs),size(beta_cond_subs,2),size(beta_cond_subs,1)*size(beta_cond_subs,3));
data_RDM = reshape(mean(conds_9,1),size(conds_9,2),1*size(conds_9,3));
DisMatrix = squareform(pdist(data_RDM, 'mahal',covdiag(data_RDM)));
% % data_RDM=smoothRDM(data_RDM,9) %*************************
% DisMatrix = squareform((pdist(data_RDM, 'euclidean')));

figure; set(gcf,'unit','pixels','position',[500 50 1200 200])
% 3-D
% [Y, stress,dispa] = mdscale(DisMatrix,3);
[coeff,Y, ~,~,~] = pca(data_RDM);
Y=zscore(Y);

subplot(1,4,1)
for p_i = 1:size(Y,1)
    s=scatter3(Y(p_i,1),Y(p_i,2),Y(p_i,3),'filled');
    s.SizeData = 300; s.MarkerFaceColor=colorlist(p_i*15,:)
    hold on
    text(Y(p_i,1),Y(p_i,2),Y(p_i,3),sprintf('%.0f %.0f %.0f',round(p_i)))
end
% axis equal
fill3([Y(1,1),Y(2,1),Y(3,1)],...
      [Y(1,2),Y(2,2),Y(3,2)],...
      [Y(1,3),Y(2,3),Y(3,3)],colorlist(90,:),'EdgeColor','none')
fill3([Y(4,1),Y(5,1),Y(6,1)],...
      [Y(4,2),Y(5,2),Y(6,2)],...
      [Y(4,3),Y(5,3),Y(6,3)],colorlist(130,:),'EdgeColor','none')
fill3([Y(7,1),Y(8,1),Y(9,1)],...
      [Y(7,2),Y(8,2),Y(9,2)],...
      [Y(7,3),Y(8,3),Y(9,3)],colorlist(160,:),'EdgeColor','none') 
grid on

plot3([Y(2,1),Y(5,1),Y(8,1)],...
      [Y(2,2),Y(5,2),Y(8,2)],...
      [Y(2,3),Y(5,3),Y(8,3)],'-','LineWidth',2)
% plot3([Y(3,1),Y(6,1),Y(9,1)],...
%       [Y(3,2),Y(6,2),Y(9,2)],...
%       [Y(3,3),Y(6,3),Y(9,3)],'-','LineWidth',2)
% plot3([Y(1,1),Y(4,1),Y(7,1)],...
%       [Y(1,2),Y(4,2),Y(7,2)],...
%       [Y(1,3),Y(4,3),Y(7,3)],'-','LineWidth',2)
grid on

% legend('1:rew-Hard','2:rew-Mid','3:rew-Easy','4:1E-Hard','5:1E-Mid','6:1E-Easy','7:2E-Hard','8:2E-Mid','9:2E-Easy')

xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% view(40,30)
hold off


% 2-D
subplot(1,4,2)
for p_i = 1:size(Y,1)
    s=scatter(Y(p_i,1),Y(p_i,2),'filled');
    s.SizeData = 300; s.MarkerFaceColor=colorlist(p_i*15,:)
    hold on
    text(Y(p_i,1),Y(p_i,2),sprintf('%.0f %.0f',round(p_i)))
end
% axis equal
grid on
% legend('Hard','Mid','Easy','0Er','1E','2E','Diff','Ners')
xlabel('PC1'); ylabel('PC2');
hold off

% 2-D
subplot(1,4,3)
for p_i = 1:size(Y,1)
    s=scatter(Y(p_i,1),Y(p_i,3),'filled');
    s.SizeData = 300; s.MarkerFaceColor=colorlist(p_i*15,:)
    hold on
    text(Y(p_i,1),Y(p_i,3),sprintf('%.0f %.0f',round(p_i)))
end
    plot(Y([2 3],1),Y([2 3],3),'--','LineWidth',2,'Color',colorlist(30,:));
    plot(Y([5 6],1),Y([5 6],3),'--','LineWidth',2,'Color',colorlist(90,:));
    plot(Y([8 9],1),Y([8 9],3),'--','LineWidth',2','Color',colorlist(130,:));
    plot(Y([1 4 7],1),Y([1 4 7],3),'-','LineWidth',3','Color',colorlist(130,:));
% axis equal
grid on
xlabel('PC1'); ylabel('PC3');
hold off

% 2-D
subplot(1,4,4)
for p_i = 1:size(Y,1)
    s=scatter(Y(p_i,2),Y(p_i,3),'filled');
    s.SizeData = 300; s.MarkerFaceColor=colorlist(p_i*15,:)
    hold on
    text(Y(p_i,2),Y(p_i,3),sprintf('%.0f %.0f',round(p_i)))
end

    plot(Y([2 1 3],2),Y([2 1 3],3),'-','LineWidth',2,'Color',colorlist(30,:));
    plot(Y([5 4 6],2),Y([5 4 6],3),'-','LineWidth',2,'Color',colorlist(90,:));
    plot(Y([8 7 9],2),Y([8 7 9],3),'-','LineWidth',2','Color',colorlist(130,:));
    %for 0 1 2 errs
    plot(Y([5 2 8],2),Y([5 2 8],3),'-','LineWidth',4','Color',colorlist(1,:));
grid on
xlabel('PC2'); ylabel('PC3');
hold off

print(['GLM2_Dev_Ners_Geome_MDS_RNN_pca0'], '-dsvg','-painters');










%% for cal mahalonobis
 function sigma=covdiag(x)

    [t,n]=size(x);
    meanx=mean(x);
    x=x-meanx(ones(t,1),:);

    sample=(1/t).*(x'*x);

    prior=diag(diag(sample));

    d=1/n*norm(sample-prior,'fro')^2;
    y=x.^2;
    r2=1/n/t^2*sum(sum(y'*y))-1/n/t*sum(sum(sample.^2));

    shrinkage=max(0,min(1,r2/d));
    sigma=shrinkage*prior+(1-shrinkage)*sample;
 end
