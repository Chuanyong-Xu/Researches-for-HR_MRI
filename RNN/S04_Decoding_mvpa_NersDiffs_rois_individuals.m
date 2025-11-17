%% 
clear; close all; clc
restoredefaultpath
imdir2 = '/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs'
 %% -----------------extract the  bold time series in mask after RDM confidence-----------------
nii_dir  = dir([imdir2,filesep, 'decoding_RDMforRNN', filesep, 'data_RNN*.csv']); % find the nii
nii_dir_sub_data9={}; nii_dir_sub_data={}; nii_dir_sub_data_max=[];

for subs = 1:length(nii_dir)

    for conds = 1:9
        whole_data = readtable([nii_dir(subs).folder filesep nii_dir(subs).name]);
        roi_data    = whole_data(whole_data.conds_9==conds,:); % rois
        nii_dir_sub_data9{subs}(:,conds)=mean(table2array(roi_data(:,13:44)))';
        nii_dir_sub_data9{subs}(:,conds)=nii_dir_sub_data9{subs}(:,conds);%./max(nii_dir_sub_data(:,conds));
    end
    %ners
    nii_dir_sub_data{subs,1}(:,1) = mean(nii_dir_sub_data9{subs}(:,[1 2 3]),2);
    nii_dir_sub_data{subs,1}(:,2) = mean(nii_dir_sub_data9{subs}(:,[4 5 6]),2);
    nii_dir_sub_data{subs,1}(:,3) = mean(nii_dir_sub_data9{subs}(:,[7 8 9]),2);
    %diffs
    nii_dir_sub_data{subs,2}(:,1) = mean(nii_dir_sub_data9{subs}(:,[1 4 7]),2);
    nii_dir_sub_data{subs,2}(:,2) = mean(nii_dir_sub_data9{subs}(:,[2 5 8]),2);
    nii_dir_sub_data{subs,2}(:,3) = mean(nii_dir_sub_data9{subs}(:,[3 6 9]),2);      
    %delete the null units
    nii_dir_sub_data{subs,1}(nii_dir_sub_data{subs,1}(:,1)==0,:)=[];%for original weights
    nii_dir_sub_data{subs,2}(nii_dir_sub_data{subs,2}(:,1)==0,:)=[];%for original weights
    % find the max
    for ind=1:size(nii_dir_sub_data{subs},1)
        nii_dir_sub_data_max{subs,1}(ind,1) = find(max(nii_dir_sub_data{subs,1}(ind,:))==(nii_dir_sub_data{subs,1}(ind,:)));
        nii_dir_sub_data_max{subs,2}(ind,1) = find(max(nii_dir_sub_data{subs,2}(ind,:))==(nii_dir_sub_data{subs,2}(ind,:)));   
    end
    
%         nii_dir_sub_data{subs} = (nii_dir_sub_data{subs})./max(nii_dir_sub_data{subs},[],2); % max value scale
    nii_dir_sub_data{subs,1} = (nii_dir_sub_data{subs,1}-min(nii_dir_sub_data{subs,1},[],2))./(max(nii_dir_sub_data{subs,1},[],2) - min(nii_dir_sub_data{subs,1},[],2)); %min-max scale
    nii_dir_sub_data{subs,2} = (nii_dir_sub_data{subs,2}-min(nii_dir_sub_data{subs,2},[],2))./(max(nii_dir_sub_data{subs,2},[],2) - min(nii_dir_sub_data{subs,2},[],2)); %min-max scale
    
    % prepare
    %ners
    % % % nii_dir_sub_data_subs0{1}(subs,:) = nanmean(nii_dir_sub_data{subs,1}(nii_dir_sub_data_max{subs,1}==1,:));
    % % % nii_dir_sub_data_subs1{1}(subs,:) = nanmean(nii_dir_sub_data{subs,1}(nii_dir_sub_data_max{subs,1}==2,:));
    % % % nii_dir_sub_data_subs2{1}(subs,:) = nanmean(nii_dir_sub_data{subs,1}(nii_dir_sub_data_max{subs,1}==3,:));
    % % % nii_dir_sub_data_subs0s{1}(subs,:) = nanstd(nii_dir_sub_data{subs,1}(nii_dir_sub_data_max{subs,1}==1,:));
    % % % nii_dir_sub_data_subs1s{1}(subs,:) = nanstd(nii_dir_sub_data{subs,1}(nii_dir_sub_data_max{subs,1}==2,:));
    % % % nii_dir_sub_data_subs2s{1}(subs,:) = nanstd(nii_dir_sub_data{subs,1}(nii_dir_sub_data_max{subs,1}==3,:));        
    % % % %diffs
    % % % nii_dir_sub_data_subs0{2}(subs,:) = nanmean(nii_dir_sub_data{subs,2}(nii_dir_sub_data_max{subs,2}==1,:));
    % % % nii_dir_sub_data_subs1{2}(subs,:) = nanmean(nii_dir_sub_data{subs,2}(nii_dir_sub_data_max{subs,2}==2,:));
    % % % nii_dir_sub_data_subs2{2}(subs,:) = nanmean(nii_dir_sub_data{subs,2}(nii_dir_sub_data_max{subs,2}==3,:));
    % % % nii_dir_sub_data_subs0s{2}(subs,:) = nanstd(nii_dir_sub_data{subs,2}(nii_dir_sub_data_max{subs,2}==1,:));
    % % % nii_dir_sub_data_subs1s{2}(subs,:) = nanstd(nii_dir_sub_data{subs,2}(nii_dir_sub_data_max{subs,2}==2,:));
    % % % nii_dir_sub_data_subs2s{2}(subs,:) = nanstd(nii_dir_sub_data{subs,2}(nii_dir_sub_data_max{subs,2}==3,:));        


    % nii_dir_sub_data_subs1{1}(subs,:) = nanmean(nii_dir_sub_data{subs,1}(nii_dir_sub_data_max{subs,1}==2 & nii_dir_sub_data_max{subs,2}==2,:)); %errors data, mu: E1, D2
    nii_dir_sub_data_subs1{1}(subs,:) = nanmean(nii_dir_sub_data{subs,1}(nii_dir_sub_data_max{subs,1}==2 ,:)); %errors data, mu: E1    
    nii_dir_sub_data_subs1s{1}(subs,:) = nanstd(nii_dir_sub_data{subs,1}(nii_dir_sub_data_max{subs,1}==2 ,:)); %sigma: E1

    nii_dir_sub_data_subs0{2}(subs,:) = nanmean(nii_dir_sub_data{subs,2}(nii_dir_sub_data_max{subs,2}==2 ,:)); %Diifculty data, mu: D2
    nii_dir_sub_data_subs0s{2}(subs,:) = nanstd(nii_dir_sub_data{subs,2}(nii_dir_sub_data_max{subs,2}==2 ,:));

end

    %% plot the gaussan
    %ners, middle
    means       = (nanmean(nii_dir_sub_data_subs1{1}));
    sigma       = (nanmean(nii_dir_sub_data_subs1s{1}));
    errs        = linspace(0,2,3);
    xFine_ners       = linspace(min(-5), max(5),100)';
    gauss       = fittype('A*exp(-((x-mu) .^2)/(2*sigma^2))','coeff',{'A','mu','sigma'});
    optsG_ners       = fitoptions('Method', 'NonlinearLeastSquares', 'Start',...
        [max(means), errs(means==max(means)),1], 'Lower',[0,1,0],'Upper',[1,1,Inf],'Display','off');
    [cfung_ners, gofg_ners] = fit(errs',means',gauss,optsG_ners); %*****************

    %diffs, middle
    means       = (nanmean(nii_dir_sub_data_subs0{2}));
    sigma       = (nanmean(nii_dir_sub_data_subs0s{2}));
    errs        = linspace(0,2,3);
    xFine_diffs       = linspace(min(-5), max(5),100)';
    gauss       = fittype('A*exp(-((x-mu) .^2)/(2*sigma^2))','coeff',{'A','mu','sigma'});
    optsG_diffs       = fitoptions('Method', 'NonlinearLeastSquares', 'Start',...
        [max(means), errs(means==max(means)),1], 'Lower',[0,1,0],'Upper',[1,1,Inf],'Display','off');
    [cfung_diffs, gofg_diffs] = fit(errs',means',gauss,optsG_diffs); %*****************


%% 3-D Gaussian
% 假设要画一个纯粹的二维正态分布曲面，横轴和纵轴都一起动。
%(select the units respond highest to the middle item)
x = linspace(min(-5), max(5),100)';      % X
y = linspace(min(-5), max(5),100)';      % Y 

[X_2d, Y_2d] = meshgrid(x, y);

% 2-D Gaussian
% Z = exp( -(((X_2d - cfung_ners.mu).^2)/(2*cfung_ners.sigma^2) + ((Y_2d - cfung_diffs.mu).^2)/(2*cfung_diffs.sigma^2)) );

% 2-D separable Gaussian with amplitudes
Rx = cfung_ners.A  * exp(-((X_2d - cfung_ners.mu).^2)/(2*cfung_ners.sigma^2));
Ry = cfung_diffs.A * exp(-((Y_2d - cfung_diffs.mu).^2)/(2*cfung_diffs.sigma^2));
Z  = Rx .* Ry;   %


addpath(genpath('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/matlab-special-heatmap'));
colorlist = flip(slanCM(21));%4, 5, 21
% colors=[colorlist(20,:); colorlist(50,:); colorlist(90,:)];
figure;
set(gcf,'unit','pixels','position',[200 50 450 400]);
surf(X_2d, Y_2d, Z);
shading interp;
colormap(colorlist(1:230,:));
colorbar
view(50,30);
hold on
% 添加网格线 - X方向
num_lines_x = 11;  % X方向的线条数量
x_lines = linspace(min(x), max(x), num_lines_x);
for i = 1:length(x_lines)
    % 固定X值，绘制Y方向的曲线
    x_val = x_lines(i);
    % Z_line = exp( -(((x_val - cfung_ners.mu).^2)/(2*cfung_ners.sigma^2) + ((y - cfung_diffs.mu).^2)/(2*cfung_diffs.sigma^2)) );
    Rx = cfung_ners.A  * exp(-((x_val - cfung_ners.mu).^2)/(2*cfung_ners.sigma^2));
    Ry = cfung_diffs.A * exp(-((y - cfung_diffs.mu).^2)/(2*cfung_diffs.sigma^2));
    Z_line  = Rx .* Ry;   %
    plot3(ones(size(y)) * x_val, y, Z_line, 'w-', 'LineWidth', .2);
end
% 添加网格线 - Y方向
num_lines_y = 11;  % Y方向的线条数量
y_lines = linspace(min(y), max(y), num_lines_y);
for i = 1:length(y_lines)
    % 固定Y值，绘制X方向的曲线
    y_val = y_lines(i);
    % Z_line = exp( -(((x - cfung_ners.mu).^2)/(2*cfung_ners.sigma^2) + ((y_val - cfung_diffs.mu).^2)/(2*cfung_diffs.sigma^2)) );
    Rx = cfung_ners.A  * exp(-((x - cfung_ners.mu).^2)/(2*cfung_ners.sigma^2));
    Ry = cfung_diffs.A * exp(-((y_val - cfung_diffs.mu).^2)/(2*cfung_diffs.sigma^2));
    Z_line  = Rx .* Ry;   %
    plot3(x, ones(size(x)) * y_val, Z_line, 'w-', 'LineWidth', .2);
end
hold off;


xlabel('Errors', 'FontSize',14);
xticks([0 1 2]); xticklabels({'E0','E1','E2'});
% xticks([0 1 2]); xticklabels({'E0','E1','E2'});

ylabel('Difficulty', 'FontSize',14);
yticks([0 1 2]); yticklabels({'D1','D2','D3'});
% xticks([0 1 2]); xticklabels({'E0','E1','E2'});

zlabel('Response', 'FontSize',14);
title('units sensitive to E1 & D2', 'FontSize',18,'FontWeight','Normal');

% print(['figures/' 'Ners_Diffs_gau3D'], '-dsvg','-painters');
saveas(gcf, ['figures/' 'Ners_Diffs_gau3D.svg']);

