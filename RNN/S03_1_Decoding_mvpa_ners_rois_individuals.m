%% 
clear; close all; clc
restoredefaultpath
imdir2 = '/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs'
 %% -----------------extract the  bold time series in mask after RDM confidence-----------------
nii_dir  = dir([imdir2,filesep, 'decoding_RDMforRNN', filesep, 'data_RNN*.csv']); % find the nii

nii_dir_sub_data9={}; nii_dir_sub_data={}; nii_dir_sub_data_max=[];

figure; set(gcf,'unit','pixels','position',[200 50 1600 1600]);

    for subs = 1:length(nii_dir)

        for conds = 1:9
            whole_data = readtable([nii_dir(subs).folder filesep nii_dir(subs).name]);
            roi_data    = whole_data(whole_data.conds_9==conds,:); % rois
            nii_dir_sub_data9{subs}(:,conds)=mean(table2array(roi_data(:,13:44)))';
            nii_dir_sub_data9{subs}(:,conds)=nii_dir_sub_data9{subs}(:,conds);%./max(nii_dir_sub_data(:,conds));
        end
        nii_dir_sub_data{subs}(:,1) = mean(nii_dir_sub_data9{subs}(:,[1 2 3]),2);
        nii_dir_sub_data{subs}(:,2) = mean(nii_dir_sub_data9{subs}(:,[4 5 6]),2);
        nii_dir_sub_data{subs}(:,3) = mean(nii_dir_sub_data9{subs}(:,[7 8 9]),2);

        nii_dir_sub_data{subs}(nii_dir_sub_data{subs}(:,1)==0,:)=[];%for original weights
        
        % find the max
        for ind=1:size(nii_dir_sub_data{subs},1)
            nii_dir_sub_data_max{subs}(ind,1) = find(max(nii_dir_sub_data{subs}(ind,:))==(nii_dir_sub_data{subs}(ind,:)));
        end
        
%         nii_dir_sub_data{subs} = (nii_dir_sub_data{subs})./max(nii_dir_sub_data{subs},[],2); % max value scale
        nii_dir_sub_data{subs} = (nii_dir_sub_data{subs}-min(nii_dir_sub_data{subs},[],2))./(max(nii_dir_sub_data{subs},[],2) - min(nii_dir_sub_data{subs},[],2)); %min-max scale
        
        % for plotting individuals
        subplot(5,6,subs); hold on
        addpath(genpath('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/matlab-special-heatmap'));
        colorlist = flip(slanCM(2));%ner: 2 dif :77
        colors=[colorlist(20,:); colorlist(50,:); colorlist(90,:)];
            mu0 = mean(nii_dir_sub_data{subs}(nii_dir_sub_data_max{subs}==1,:),1);
            plot([0:2],mu0,'-o','LineWidth',2,'Color',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:));    
            mu1 = mean(nii_dir_sub_data{subs}(nii_dir_sub_data_max{subs}==2,:),1);
            plot([0:2],mu1,'-o','LineWidth',2,'Color',colors(2,:),'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:));    
            mu2 = mean(nii_dir_sub_data{subs}(nii_dir_sub_data_max{subs}==3,:),1);
            plot([0:2],mu2,'-o','LineWidth',2,'Color',colors(3,:),'MarkerFaceColor',colors(3,:),'MarkerEdgeColor',colors(3,:));

        legend({'0er','1er','2er'}, 'Location','NorthEast'); 
        xlabel('Number of errors'); ylabel('Weights');
        hold off
        
        % prepare
        nii_dir_sub_data_subs0(subs,:) = nanmean(nii_dir_sub_data{subs}(nii_dir_sub_data_max{subs}==1,:));
        nii_dir_sub_data_subs1(subs,:) = nanmean(nii_dir_sub_data{subs}(nii_dir_sub_data_max{subs}==2,:));
        nii_dir_sub_data_subs2(subs,:) = nanmean(nii_dir_sub_data{subs}(nii_dir_sub_data_max{subs}==3,:));
        
        nii_dir_sub_data_subs0s(subs,:) = nanstd(nii_dir_sub_data{subs}(nii_dir_sub_data_max{subs}==1,:));
        nii_dir_sub_data_subs1s(subs,:) = nanstd(nii_dir_sub_data{subs}(nii_dir_sub_data_max{subs}==2,:));
        nii_dir_sub_data_subs2s(subs,:) = nanstd(nii_dir_sub_data{subs}(nii_dir_sub_data_max{subs}==3,:));        
    end
print(['figures/' 'Ners_weights_density_rois_subs_' num2str(subs)],'-dsvg','-painters');
% close all
        %% plot original normalized weights
        addpath(genpath('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/matlab-special-heatmap'));
        colorlist = flip(slanCM(2));%ner: 2 dif :77
        colors=[colorlist(20,:); colorlist(50,:); colorlist(90,:)];

    figure; hold on
    set(gcf,'unit','pixels','position',[200 50 400 200]);
            mu0 = nanmean(nii_dir_sub_data_subs0);
            plot([0:2],mu0,'-o','LineWidth',2,'Color',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:)); 
            mu1 = nanmean(nii_dir_sub_data_subs1);
            plot([0:2],mu1,'-o','LineWidth',2,'Color',colors(2,:),'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:))    
            mu2 = nanmean(nii_dir_sub_data_subs2);
            plot([0:2],mu2,'-o','LineWidth',2,'Color',colors(3,:),'MarkerFaceColor',colors(3,:),'MarkerEdgeColor',colors(3,:))  

        legend({'0er','1er','2er'}, 'Location','NorthEast'); 
        xlabel('Number of errors'); ylabel('Normalized Weights');
        xticks([0 1 2]); xticklabels({'0er','1er','2er'});
        hold off
print(['figures/' 'Ners_weights_density_rois_subsMean'],'-dsvg','-painters');
close all 
    %% plot the gaussan
    figure; hold on
    set(gcf,'unit','pixels','position',[200 50 400 200]);
    set(0, 'DefaultAxesFontName', 'Helvetica');   % Fonttype for axis
    set(0, 'DefaultTextFontName', 'Helvetica');   % Fonttype for text
    set(gca,'FontSize',14,'tickdir','out')
    
    means       = (nanmean(nii_dir_sub_data_subs0));
    sigma       = (nanmean(nii_dir_sub_data_subs0s));
%     semw        = sigma./sqrt(length(nii_dir_sub_data(nii_dir_sub_data_max==1,:)));
    errs        = linspace(0,2,3);
    xFine       = linspace(min(-5), max(5),200000)';
    gauss       = fittype('A*exp(-((x-mu) .^2)/(2*sigma^2))','coeff',{'A','mu','sigma'});
    optsG       = fitoptions('Method', 'NonlinearLeastSquares', 'Start',...
        [max(means), errs(means==max(means)),1], 'Lower',[0,0,0],'Upper',[1,0,Inf],'Display','off');%,'Weights',1./(semw.^2)
    [cfung, gofg] = fit(errs',means',gauss,optsG);%*****************
    ygauss      = cfung.A * exp(-((xFine-cfung.mu) .^2)/(2*cfung.sigma^2));
    plot(xFine,ygauss,'LineWidth',2,'Color',colors(1,:)); hold on
    plot(errs, means,'o','MarkerSize',6,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:))
%  
    means       = (nanmean(nii_dir_sub_data_subs1));
    sigma       = (nanmean(nii_dir_sub_data_subs1s));
    errs        = linspace(0,2,3);
    xFine       = linspace(min(-5), max(5),200000)';
    gauss       = fittype('A*exp(-((x-mu) .^2)/(2*sigma^2))','coeff',{'A','mu','sigma'});
    optsG       = fitoptions('Method', 'NonlinearLeastSquares', 'Start',...
        [max(means), errs(means==max(means)),1], 'Lower',[0,0,0],'Upper',[1,1,Inf],'Display','off');
    [cfung, gofg] = fit(errs',means',gauss,optsG); %*****************
    ygauss      = cfung.A * exp(-((xFine-cfung.mu) .^2)/(2*cfung.sigma^2));
    plot(xFine,ygauss,'LineWidth',2,'Color',colors(2,:)); hold on
    plot(errs, means,'o','MarkerSize',6,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:))
%    
    means       = (nanmean(nii_dir_sub_data_subs2));
    sigma       = (nanmean(nii_dir_sub_data_subs2s));
    errs        = linspace(0,2,3);
    xFine       = linspace(min(-5), max(5),200000)';
    gauss       = fittype('A*exp(-((x-mu) .^2)/(2*sigma^2))','coeff',{'A','mu','sigma'});
    optsG       = fitoptions('Method', 'NonlinearLeastSquares', 'Start',...
        [max(means), errs(means==max(means)),1], 'Lower',[0,0,0],'Upper',[1,2,Inf],'Display','off');
    [cfung, gofg] = fit(errs',means',gauss,optsG); %*****************
    ygauss      = cfung.A * exp(-((xFine-cfung.mu) .^2)/(2*cfung.sigma^2));
    plot(xFine,ygauss,'LineWidth',2,'Color',colors(3,:)); hold on
    plot(errs, means,'o','MarkerSize',6,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor',colors(3,:))
    
    legend({'0er','0er','1er','1er','2er','2er'}, 'Location','NorthWest'); 
    xlabel('Number of errors'); ylabel('Normalized Weights')
    xticks([0 1 2]); xticklabels({'0er','1er','2er'});
    hold off
print(['figures/' 'Ners_weights_density_rois_subsMean_gau'],'-dsvg','-painters');
% close all


