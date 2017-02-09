%% Determine linear transformation between two datasets based on fiducial trajectory

% Input : 2D localization data

% Output: dx and dy to correct the A750 channel

clear, clc, close all
%% Load the two datasets 

%  Load the cleaned version with beads

% Format 1: bstore 2: TS

formatC1 = 2;
formatC2 = 1;

path_Ch1            = 'Z:\Christian-Sieben\data_HTP\2016-09-22_A549_EGFR_SNA\locResults\A549_SNA_A647_18';
filename_locs_Ch1   = 'A549_SNA_A647_18_MMStack_Pos0_locResults_TS_DC'; 

path_Ch2            = 'Z:\Christian-Sieben\data_HTP\2016-09-22_A549_EGFR_SNA\locResults\A549_EGFR_A750_19';
filename_locs_Ch2   = 'A549_EGFR_A750_19_MMStack_Pos0_locResults_cleaned'; 

cd(path_Ch1);
locs_Ch1=dlmread([filename_locs_Ch1 '.dat'],',',1,0);
cd(path_Ch2);
locs_Ch2=dlmread([filename_locs_Ch2 '.dat'],',',1,0);

% Load Header Data

cd(path_Ch1);
file        = fopen([filename_locs_Ch1 '.dat']);
line        = fgetl(file);
h           = regexp( line, ',', 'split' );

if formatC1 == 1;

xCol        = strmatch('x [nm]',h);
yCol        = strmatch('y [nm]',h);
frameCol    = strmatch('frame',h);
photonsCol  = strmatch('intensity [photon]',h);

else 
   

xCol        = strmatch('"x [nm]"',h);
yCol        = strmatch('"y [nm]"',h);
frameCol    = strmatch('"frame"',h);
photonsCol  = strmatch('"intensity [photon]"',h);

end

cd(path_Ch2);
file        = fopen([filename_locs_Ch2 '.dat']);
line        = fgetl(file);
h           = regexp( line, ',', 'split' );


if formatC2 == 1;   
    
xCol2        = strmatch('x [nm]',h);
yCol2        = strmatch('y [nm]',h);
dxCol2       = strmatch('dx',h);
dyCol2       = strmatch('dy',h);
frameCol2    = strmatch('frame',h);
photonsCol2  = strmatch('intensity [photon]',h);

else 

xCol2        = strmatch('"x [nm]"',h);
yCol2        = strmatch('"y [nm]"',h);
frameCol2    = strmatch('"frame"',h);
photonsCol2  = strmatch('"intensity [photon]"',h);

end
    
fprintf('\n -- Data loaded --\n')

%% Plot an overlay of a subset 

minFrame = 1e3;
maxFrame = 5e4;

figure
scatter(locs_Ch1(minFrame:maxFrame,xCol),locs_Ch1(minFrame:maxFrame,yCol),1,'green'); hold on;
scatter(locs_Ch2(minFrame:maxFrame,xCol),locs_Ch2(minFrame:maxFrame,yCol),1,'red');

%% Plot only one channel  

figure
scatter(locs_Ch1(minFrame:maxFrame,xCol),locs_Ch1(minFrame:maxFrame,yCol),1);


%% Plot the data, Show 2D histogram, Select Au Fuducial using rectangular selection

close all
pxlsize = 500; 

heigth=round((max(locs_Ch1(:,yCol))-min(locs_Ch1(:,yCol)))/pxlsize);
width=round((max(locs_Ch1(:,xCol))-min(locs_Ch1(:,xCol)))/pxlsize);

figure('Position',[650 400 500 500])
im=hist3([locs_Ch2(:,xCol),locs_Ch2(:,yCol)],[width heigth]); % heigth x width
imagesc(imrotate(im,90),[10 10000]);
colormap('hot');
colorbar
rect = getrect; % rect = [xmin ymin width height];
close all;

fprintf('\n -- ROI selected --\n')

% Select ROI with Au Fiducial

xmin = min(locs_Ch1(:,xCol))+ rect(:,1) * pxlsize;
ymin = max(locs_Ch1(:,yCol)) - rect(:,2)* pxlsize - (rect(:,4)*pxlsize) ;
xmax = xmin + (rect(:,3)* pxlsize);
ymax = ymin + rect(:,4) * pxlsize;

% Select ROI for Channel 1

vx=find(locs_Ch1(:,xCol)>xmin & locs_Ch1(:,xCol)<xmax);
subset1=locs_Ch1(vx,1:end);
vy=find(subset1(:,yCol)>ymin & subset1(:,yCol)<ymax);
subset2=subset1(vy,1:end);

figure('Position',[100 400 500 500],'Name','Drift Ch1')
scatter(subset2(:,xCol),subset2(:,yCol),1)

% Select ROI for Channel 2

vx=find(locs_Ch2(:,xCol)>xmin & locs_Ch2(:,xCol)<xmax);
subset3=locs_Ch2(vx,1:end);
vy=find(subset3(:,yCol)>ymin & subset3(:,yCol)<ymax);
subset4=subset3(vy,1:end);

figure('Position',[1200 400 500 500],'Name','Drift Ch2')
scatter(subset4(:,xCol),subset4(:,yCol),1)

fprintf('\n -- Plotted selected ROI  --\n')

%% Refine ROI for Ch1

xmin = 3.8e4;
xmax = 3.83e4;
ymin = 2.66e4;
ymax = 2.69e4;

% xmin=min(locs(:,3));
% xmax=max(locs(:,3)); 
% ymin=min(locs(:,4));
% ymax=max(locs(:,4));

vx=find(locs_Ch1(:,xCol)>xmin & locs_Ch1(:,xCol)<xmax);
subset1=locs_Ch1(vx,1:end);
vy=find(subset1(:,yCol)>ymin & subset1(:,yCol)<ymax);
subset2=subset1(vy,1:end);

figure('Position',[1200 400 500 500])
scatter(subset2(:,xCol),subset2(:,yCol),1)

fprintf('\n -- Plotted selected ROI  --\n')

%% Refine ROI for Ch2

xmin = 3.8e4;
xmax = 3.83e4;
ymin = 2.66e4;
ymax = 2.69e4;


% xmin=min(locs(:,3));
% xmax=max(locs(:,3)); 
% ymin=min(locs(:,4));
% ymax=max(locs(:,4));

vx=find(locs_Ch2(:,xCol)>xmin & locs_Ch2(:,xCol)<xmax);
subset3=locs_Ch2(vx,1:end);
vy=find(subset3(:,yCol)>ymin & subset3(:,yCol)<ymax);
subset4=subset3(vy,1:end);

figure('Position',[1200 400 500 500])
scatter(subset4(:,xCol),subset4(:,yCol),1)

fprintf('\n -- Plotted selected ROI  --\n')


%% Plot XY Scatter and XY drift

close all

figure
scatter(subset2(:,xCol),subset2(:,yCol),2,'black');hold on;
scatter(subset4(:,xCol),subset4(:,yCol),2,'green');
legend('Channel 1','Channel 2')

figure('Position',[200 200 600 600],'Name','Drift Ch1 and Ch2')
h=gcf;
set(h,'PaperOrientation','landscape');

subplot(2,2,1)
scatter(subset2(:,frameCol),(subset2(:,xCol)-min(subset2(:,xCol))),1);
title('X drift Ch1')
xlabel('time [frame]')
ylabel('x position [nm]')
box on;

subplot(2,2,2)
scatter(subset2(:,frameCol),(subset2(:,yCol)-min(subset2(:,yCol))),1);
title('Y drift Ch1')
xlabel('time [frame]')
ylabel('y position [nm]')
box on;

subplot(2,2,3)
scatter(subset4(:,frameCol),(subset4(:,xCol)-min(subset4(:,xCol))),1);
title('X drift Ch2')
xlabel('time [frame]')
ylabel('x position [nm]')
box on;

subplot(2,2,4)
scatter(subset4(:,frameCol),(subset4(:,yCol)-min(subset4(:,yCol))),1);
title('Y drift Ch2')
xlabel('time [frame]')
ylabel('y position [nm]')
box on;

%% Fit the drift trajectory

close all;

%%%%%%%

minFrame  = 500; % this are the locs
maxFrame  = 1000;

minFrames = 1;
maxFrames = 15000;

%%%%%%% Channel 1 %%%%%%% 

xvalues = zeros(maxFrames,1);
xvalues(subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,frameCol)) = subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,frameCol);
xvalues(1:minFrames) = 1:minFrames;

yvalues = zeros(maxFrames,1);
yvalues(subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,frameCol)) = (subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,xCol)-(sum(subset2(minFrame:maxFrame,xCol))/length(subset2(minFrame:maxFrame,xCol))));
yvalues(1:minFrames) = mean(yvalues(minFrames:minFrames+100));

[fxC1,gofxC1,outputxC1] = fit(xvalues,yvalues,'poly8');
figure
scatter(xvalues, yvalues);hold on;
plot(fxC1)

xvalues = zeros(maxFrames,1);
xvalues(subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,frameCol)) = subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,frameCol);
xvalues(1:minFrames) = 1:minFrames;

yvalues = zeros(maxFrames,1);
yvalues(subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,frameCol)) = (subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,yCol)-(sum(subset2(minFrame:maxFrame,yCol))/length(subset2(minFrame:maxFrame,yCol))));
yvalues(1:minFrames) = mean(yvalues(minFrames:minFrames+100));

[fyC1,gofyC1,outputyC1] = fit(xvalues,yvalues,'poly8');
figure
scatter(xvalues, yvalues);hold on;
plot(fyC1)

%%%%%%% Channel 2 %%%%%%% 

xvalues = zeros(maxFrames,1);
xvalues(subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,frameCol)) = subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,frameCol);
xvalues(1:minFrames) = 1:minFrames;

yvalues = zeros(maxFrames,1);
yvalues(subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,frameCol)) = (subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,xCol)-(sum(subset4(minFrame:maxFrame,xCol))/length(subset4(minFrame:maxFrame,xCol))));
yvalues(1:minFrames) = mean(yvalues(minFrames:minFrames+100));

[fxC2,gofxC2,outputxC2] = fit(xvalues,yvalues,'poly8');
figure
scatter(xvalues, yvalues);hold on;
plot(fxC2)

xvalues = zeros(maxFrames,1);
xvalues(subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,frameCol)) = subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,frameCol);
xvalues(1:minFrames) = 1:minFrames;

yvalues = zeros(maxFrames,1);
yvalues(subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,frameCol)) = (subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,yCol)-(sum(subset4(minFrame:maxFrame,yCol))/length(subset4(minFrame:maxFrame,yCol))));
yvalues(1:minFrames) = mean(yvalues(minFrames:minFrames+100));

[fyC2,gofyC2,outputyC2] = fit(xvalues,yvalues,'poly8');
figure
scatter(xvalues,yvalues);hold on;
plot(fyC2)



% [fxC1,gofxC1,outputxC1] = fit(subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,frameCol),(subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,xCol)-(sum(subset2(minFrame:maxFrame,xCol))/length(subset2(minFrame:maxFrame,xCol)))),'poly8'); % fit in nm, fx
% [fyC1,gofyC1,outputyC1] = fit(subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,frameCol),(subset2(subset2(:,frameCol)>minFrames & subset2(:,frameCol)<maxFrames,yCol)-(sum(subset2(minFrame:maxFrame,yCol))/length(subset2(minFrame:maxFrame,yCol)))),'poly8'); % fit in nm, fy
% 
% [fxC2,gofxC2,outputxC2] = fit(subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,frameCol),(subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,xCol)-(sum(subset4(minFrame:maxFrame,xCol))/length(subset4(minFrame:maxFrame,xCol)))),'poly8'); % fit in nm, fx
% [fyC2,gofyC2,outputyC2] = fit(subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,frameCol),(subset4(subset4(:,frameCol)>minFrames & subset4(:,frameCol)<maxFrames,yCol)-(sum(subset4(minFrame:maxFrame,yCol))/length(subset4(minFrame:maxFrame,yCol)))),'poly8'); % fit in nm, fy
% 

% [fxC1,gofxC1,outputxC1] = fit(subset2(:,frameCol),(subset2(:,xCol)-(sum(subset2(minFrame:maxFrame,xCol))/length(subset2(minFrame:maxFrame,xCol)))),'poly8'); % fit in nm, fx
% [fyC1,gofyC1,outputyC1] = fit(subset2(:,frameCol),(subset2(:,yCol)-(sum(subset2(minFrame:maxFrame,yCol))/length(subset2(minFrame:maxFrame,yCol)))),'poly8'); % fit in nm, fy

% [fxC2,gofxC2,outputxC2] = fit(subset4(:,frameCol),(subset4(:,xCol)-min(subset4(:,xCol))),'poly8'); % fit in nm, fx
% [fyC2,gofyC2,outputyC2] = fit(subset4(:,frameCol),(subset4(:,yCol)-min(subset4(:,yCol))),'poly8'); % fit in nm, fy

% [fxC2,gofxC2,outputxC2] = fit(subset4(:,frameCol),(subset4(:,xCol)-(sum(subset4(minFrame:maxFrame,xCol))/length(subset4(minFrame:maxFrame,xCol)))),'poly8'); % fit in nm, fx
% [fyC2,gofyC2,outputyC2] = fit(subset4(:,frameCol),(subset4(:,yCol)-(sum(subset4(minFrame:maxFrame,yCol))/length(subset4(minFrame:maxFrame,yCol)))),'poly8'); % fit in nm, fy


figure('Position',[900 200 600 600],'Name','Drift fitting both Channels')
h=gcf;
set(h,'PaperOrientation','portrait');

subplot(2,2,1)
scatter(subset2(:,frameCol),(subset2(:,xCol)-(sum(subset2(minFrame:maxFrame,xCol))/length(subset2(minFrame:maxFrame,xCol)))),1);hold on;
h = plot(fxC1);set(h,'LineWidth',2);
axis([0 max(locs_Ch1(:,frameCol)) -500 500])
title('X drift Ch1')
xlabel('time [frames]')
ylabel('x position [nm]')
legend('data','fitted curve','Location','northeast');
legend('boxon');
box on;

subplot(2,2,2)
scatter(subset2(:,frameCol),(subset2(:,yCol)-(sum(subset2(minFrame:maxFrame,yCol))/length(subset2(minFrame:maxFrame,yCol)))),1);hold on;
h = plot(fyC1);set(h,'LineWidth',2);
axis([0 max(locs_Ch1(:,frameCol)) -500 500])
title('Y drift Ch1')
xlabel('time [frames]')
ylabel('y position [nm]')
legend('data','fitted curve','Location','northeast');
legend('boxon');
box on;


subplot(2,2,3)
scatter(subset4(:,frameCol),(subset4(:,xCol)-(sum(subset4(minFrame:maxFrame,xCol))/length(subset4(minFrame:maxFrame,xCol)))),1);hold on;
h = plot(fxC2);set(h,'LineWidth',2);
axis([0 max(locs_Ch2(:,frameCol)) -500 500])
title('X drift Ch2')
xlabel('time [frames]')
ylabel('x position [nm]')
legend('data','fitted curve','Location','northeast');
legend('boxon');
box on;


subplot(2,2,4)
scatter(subset4(:,frameCol),(subset4(:,yCol)-(sum(subset4(minFrame:maxFrame,yCol))/length(subset4(minFrame:maxFrame,yCol)))),1);hold on;
h = plot(fyC2);set(h,'LineWidth',2);
axis([0 max(locs_Ch2(:,frameCol)) -500 500])
title('Y drift Ch2')
xlabel('time [frames]')
ylabel('y position [nm]')
legend('data','fitted curve','Location','northeast');
legend('boxon');
box on;


% Calculate RMSE

pdx = fitdist(outputxC1.residuals,'normal') % Distribution of residuals
pdy = fitdist(outputyC1.residuals,'normal') % Distribution of residuals

%% Calculate deviation for each frame

close all;

tic

dx=[]; % this value must be sustracted from Ch2s xcoordinate
dy=[]; % this value must be sustracted from Ch2s ycoordinate


offx = (sum(subset2(minFrame:maxFrame,xCol))/length(subset2(minFrame:maxFrame,xCol)))-(sum(subset4(minFrame:maxFrame,xCol))/length(subset4(minFrame:maxFrame,xCol)));
offy = (sum(subset2(minFrame:maxFrame,yCol))/length(subset2(minFrame:maxFrame,yCol)))-(sum(subset4(minFrame:maxFrame,yCol))/length(subset4(minFrame:maxFrame,yCol)));

for frame = minFrames:maxFrames;%(locs_Ch2(:,frameCol));

dx(frame,1) = fxC1(frame) - fxC2(frame) - offx ; % in nm
dy(frame,1) = fyC1(frame) - fyC2(frame) - offy ; % in nm

end

toc


%%  Test Correction on bead image

subset4C=[];

j=1;

for i = 1:length(subset4);
       
frame = subset4(i,frameCol);

if frame > 8000 | frame < minFrames;
    
    subset4C(j,1) = (subset4(i,xCol));
    subset4C(j,2) = (subset4(i,yCol));
    
else
    
subset4C(j,1) = (subset4(i,xCol)) - (dx(frame));
subset4C(j,2) = (subset4(i,yCol)) - (dy(frame));

% subset4C(j,1)=(subset2(i,xCol))+abs(dx(frame));
% subset4C(j,2)=(subset2(i,yCol))-abs(dy(frame));
end
clear frame

j=j+1;

end

figure
scatter(subset2(:,xCol),subset2(:,yCol),2,'black');hold on;
scatter(subset4C(:,1),subset4C(:,2),2,'red');
scatter(subset4(:,xCol),subset4(:,yCol),2,'green');
legend('Ch1','Ch2 corrected','Ch2 ');
box on;


%% 
cd(path_Ch2);
save('devX.mat','dx');
save('devY.mat','dy');
minFrames
maxFrames



%%%% Apply Correction on full Ch2 dataset
% 
% filename_uncorr_Ch2='2016-08-19_MitoRNAGran_A750_FOV_2_MMStack_Pos0_locResults_DC'; % filename of TS output file
% uncorr_Ch2=dlmread([filename_uncorr_Ch2 '.dat'],',',1,0);
% 
% corr_Ch2 = uncorr_Ch2;
% 
% %% 
% 
% 
% for i=1:length(corr_Ch2(:,frameCol));
%     
% frame=corr_Ch2(i,frameCol);
% 
% corr_Ch2(i,xCol)=(corr_Ch2(i,xCol))-abs(dx(frame)-50);
% corr_Ch2(i,yCol)=(corr_Ch2(i,yCol))+abs(dy(frame)+50);
% 
% clear frame
% 
% end
% 
% dlmwrite([filename_uncorr_Ch2 '_corr.dat'],corr_Ch2);
% 




