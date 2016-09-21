% The script takes input from LAMA CBC analysis 
% and gives some visualization options

% CBC analysis performed in Lama (v.1606), Single Molecule Biophysics, University of Frankfurt

clear, clc, close all

%% Load the CBC results 

name_locs_Ch1 = 'origin_cbc_roi'; 
name_locs_Ch2 = 'partner_cbc_roi';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

locsC1=dlmread([name_locs_Ch1 '.txt'],'',2,0);
locsC2=dlmread([name_locs_Ch2 '.txt'],'',2,0);

xCol        = 1;
yCol        = 2;
frameCol    = 3;
CBC         = 4;

%% Figure option #1: PLot overview figure Localizations and CBC

figure('Position',[100 100 1000 900])

subplot(3,2,1);
scatter(locsC1(:,xCol),locsC1(:,yCol),1,'red');hold on;
title('Locs Channel 1');
xlabel('x [nm]');
ylabel('y [nm]');
box on;
axis square;

subplot(3,2,2);
scatter(locsC2(:,xCol),locsC2(:,yCol),1,'green');hold on;
title('Locs Channel 2');
xlabel('x [nm]');
ylabel('y [nm]');
box on;
axis square;

subplot(3,2,3);
scatter(locsC1(:,xCol),locsC1(:,yCol),5,locsC1(:,CBC),'filled');hold on;
title('CBC of Channel 1');
xlabel('x [nm]');
ylabel('y [nm]');
box on;
axis square;
colormap jet; colorbar

subplot(3,2,4);
scatter(locsC2(:,xCol),locsC2(:,yCol),5,locsC2(:,CBC),'filled');hold on;
title('CBC of Channel 2');
xlabel('x [nm]');
ylabel('y [nm]');
box on;
axis square;
colormap jet; colorbar

subplot(3,2,5);
scatter(locsC1((locsC1(:,CBC)>=0.7),xCol),locsC1((locsC1(:,CBC)>=0.7),yCol),1,'red');hold on;
title('CBC Ch1 > 0.7');
xlabel('x [nm]');
ylabel('y [nm]');
% legend('CBC Ch1 > 0.7','Ch2');
box on;
axis square;

subplot(3,2,6);
scatter(locsC2((locsC2(:,CBC)>=0.7),xCol),locsC2((locsC2(:,CBC)>=0.7),yCol),1,'green');hold on;
title('CBC Ch2 > 0.7');
xlabel('x [nm]');
ylabel('y [nm]');
% legend('CBC Ch1 > 0.7','Ch2');
box on;
axis square;

%% Figure option #2: Show Overlay of both channels and CA histograms

figure('Position',[100 600 900 300])

subplot(1,3,1);
scatter(locsC1(:,xCol),locsC1(:,yCol),1,'red');hold on;
scatter(locsC2(:,xCol),locsC2(:,yCol),1,'green');hold on;
title('Color overlay');
xlabel('x [nm]');
ylabel('y [nm]');
legend('Ch1','Ch2');
box on;
axis square;


subplot(1,3,2);
bins = -1:0.1:1;
b = bar(bins, hist(locsC1(:,CBC),bins)/sum(hist(locsC1(:,CBC),bins)));
b.FaceColor = [0 0 0];
b.EdgeColor = [0.5 0.5 0.5];
b.LineWidth = 0.1;
axis([-1.1 1.1 0 0.1])
title(['C_A Histogram Channel 1 - ' num2str(sum(locsC1(:,CBC)>=0.7)/length(locsC1),'%.2f') ' colocalized']);
xlabel('C_A');
ylabel('counts');
box on;
axis square;

subplot(1,3,3);
bins = -1:0.1:1;
b = bar(bins, hist(locsC2(:,CBC),bins)/sum(hist(locsC2(:,CBC),bins)));
b.FaceColor = [0 0 0];
b.EdgeColor = [0.5 0.5 0.5];
b.LineWidth = 0.1;
axis([-1.1 1.1 0 0.1])
title(['C_A Histogram Channel 2 - ' num2str(sum(locsC2(:,CBC)>=0.7)/length(locsC2),'%.2f') ' colocalized']);
xlabel('C_A');
ylabel('counts');
box on;
axis square;

%% Figure option #3 : Calculate localization density
tic

NNdist=100; % neighbor distance in nm

idx1 = rangesearch(locsC1(:,xCol:yCol),locsC1(:,xCol:yCol),NNdist);
idx2 = rangesearch(locsC2(:,xCol:yCol),locsC2(:,xCol:yCol),NNdist);
        
for i=1:length(idx1);
    
NoNindC1(i,1)=length(idx1{i,1});     % count the total number of neighbors for each point in dataset

end

for i=1:length(idx2);
    
NoNindC2(i,1)=length(idx2{i,1});     % count the total number of neighbors for each point in dataset

end

toc
%% Figure option #3 :  Plot Density
figure('Position',[100 100 1600 800])

densityFilter = 20;

target1 = find(NoNindC1>densityFilter);
target2 = find(NoNindC2>densityFilter);
        
        ax1 = subplot(2,3,1)

scatter(locsC1(target1,1),locsC1(target1,2),4,NoNindC1(target1),'o','filled');hold on;
colormap(ax1,hot)
title(['Density Map Ch1, r = ', num2str(NNdist)]);
axis([min(locsC1(target1,1)) max(locsC1(target1,1)) min(locsC1(target1,2)) max(locsC1(target1,2))])
xlabel('nm');
ylabel('nm');
box on;


        ax2 = subplot(2,3,2)

scatter(locsC2(target2,1),locsC2(target2,2),6,NoNindC2(target2),'o','filled');hold on;
colormap(ax2,hot)
title(['Density Map Ch2, r =', num2str(NNdist)]);
axis([min(locsC2(target2,1)) max(locsC2(target2,1)) min(locsC2(target2,2)) max(locsC2(target2,2))])
xlabel('nm');
ylabel('nm');
box on;

        ax3 = subplot(2,3,4)

scatter(locsC1(target1,xCol),locsC1(target1,yCol),5,locsC1(target1,CBC),'filled');hold on;
axis([min(locsC1(target1,1)) max(locsC1(target1,1)) min(locsC1(target1,2)) max(locsC1(target1,2))])
title('CBC of Ch1');
xlabel('x [nm]');
ylabel('y [nm]');
box on;
axis square;
colorbar; colormap(ax3,jet); 

        ax4 = subplot(2,3,5)

scatter(locsC2(target2,xCol),locsC2(target2,yCol),5,locsC2(target2,CBC),'filled');hold on;
axis([min(locsC2(target2,1)) max(locsC2(target2,1)) min(locsC2(target2,2)) max(locsC2(target2,2))])
title('CBC of Ch2');
xlabel('x [nm]');
ylabel('y [nm]');
box on;
axis square;
colorbar; colormap(ax4,jet);

        subplot(2,3,3)

bins = -1:0.05:1;
b = bar(bins, hist(locsC1(target1,CBC),bins)/sum(hist(locsC1(target1,CBC),bins)));
b.FaceColor = [0 0 0];
b.EdgeColor = [0.5 0.5 0.5];
b.LineWidth = 0.1;
axis([-1.1 1.1 0 0.1])
title(['C_A score Ch1 - ' num2str(sum(locsC1(:,CBC)>=0.3)/length(locsC1),'%.2f') ' colocalized']);
xlabel('C_A');
ylabel('counts');
box on;
axis square;

        subplot(2,3,6)

bins = -1:0.05:1;
b = bar(bins, hist(locsC2(target2,CBC),bins)/sum(hist(locsC2(target2,CBC),bins)));
b.FaceColor = [0 0 0];
b.EdgeColor = [0.5 0.5 0.5];
b.LineWidth = 0.1;
axis([-1.1 1.1 0 0.1])
title(['C_A score Ch2 - ' num2str(sum(locsC2(:,CBC)>=0.3)/length(locsC2),'%.2f') ' colocalized']);
xlabel('C_A');
ylabel('counts');
box on;
axis square;





