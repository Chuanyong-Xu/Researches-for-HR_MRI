    addpath(genpath('/opt/matlab-special-heatmap'));
    colorlist = (slanCM(98));

mu = 0; sigma = 1;
    x = linspace(mu-3*sigma, mu+3*sigma, 2000);
    y = normpdf(x, mu, sigma);

    figure; set(gcf,'position',[200 50 150 50])
    set(gca, 'Color','none'); 
    hold on
    plot(x, y+0.05, 'k-', 'LineWidth', 1.5); %
    set(gca, 'XTick',[],'LineWidth',1.5);
    ax=gca; ax.YAxis.Visible='off'
    
    %
    L = 1; U = 4;
    mask = (x >= L & x <= U);
    xp = [x(mask), fliplr(x(mask))];
    yp = [zeros(1,nnz(mask)), fliplr(y(mask))+0.05];
    patch(xp, yp, colorlist(20,:), 'EdgeColor', 'none', 'FaceAlpha', .8);
    
    % 
    box off; 
%     set(gca,'TickDir','out','FontName','Helvetica','LineWidth',1.2);
%     xlabel('x'); ylabel('Density'); 

    print(gcf,'norm_fill','-dsvg','-painters');

