%% Use before CBC analysis 

% What this does:

% Visualize the FOV 
% choose ROI for CBC analysis
% export overlay images
% create a control dataset

clear, clc, close all
%% Load the two datasets 

path_Ch1            = 'Z:\Christian-Sieben\data_HTP\2016-08-19_Nucleoid_MitoRNAGran\locResults\2016-08-19_Nucleoid_A647_FOV_4';
filename_locs_Ch1   = '2016-08-19_Nucleoid_A647_FOV_4_MMStack_Pos0_locResults_cleaned'; 

path_Ch2            = 'Z:\Christian-Sieben\data_HTP\2016-08-19_Nucleoid_MitoRNAGran\locResults\2016-08-19_MitoRNAGran_A750_FOV_4';
filename_locs_Ch2   = '2016-08-19_MitoRNAGran_A750_FOV_4_MMStack_Pos0_locResults_cleaned'; 

cd(path_Ch1);
locs_Ch1=dlmread([filename_locs_Ch1 '.dat'],',',1,0);
cd(path_Ch2);
locs_Ch2=dlmread([filename_locs_Ch2 '.dat'],',',1,0);

cd(path_Ch1);
file        = fopen([filename_locs_Ch1 '.dat']);
line        = fgetl(file);
h           = regexp( line, ',', 'split' );

xCol        = strmatch('x [nm]',h);
yCol        = strmatch('y [nm]',h);
dxCol       = strmatch('dx',h);
dyCol       = strmatch('dy',h);
frameCol    = strmatch('frame',h);
photonsCol  = strmatch('intensity [photon]',h);

fprintf('\n -- Data loaded --\n')

%% Show 2 color dataset

minFrame = 1e3;
maxFrame = 1e4;

figure; hold all;
scatter(locs_Ch1(minFrame:maxFrame,xCol),locs_Ch1(minFrame:maxFrame,yCol),1,'green');
scatter(locs_Ch2(minFrame:maxFrame,xCol),locs_Ch2(minFrame:maxFrame,yCol),1,'red');

%% Select ROI

xmin = 0.5e4;
xmax = 3.5e4;
ymin = 1e4;
ymax = 2.5e4;

v1 = find(locs_Ch1(:,xCol) > xmin & locs_Ch1(:,xCol) < xmax & locs_Ch1(:,yCol) > ymin & locs_Ch1(:,yCol) < ymax);
v2 = find(locs_Ch2(:,xCol) > xmin & locs_Ch2(:,xCol) < xmax & locs_Ch2(:,yCol) > ymin & locs_Ch2(:,yCol) < ymax);

subset1=locs_Ch1(v1,1:end);
subset2=locs_Ch2(v2,1:end);

figure('Position',[100 400 500 500],'Name','Drift Ch1')
scatter(subset1(minFrame:maxFrame,xCol),subset1(minFrame:maxFrame,yCol),1,'bo')

figure('Position',[1200 400 500 500],'Name','Drift Ch2')
scatter(subset2(minFrame:maxFrame,xCol),subset2(minFrame:maxFrame,yCol),1,'bo')

fprintf('\n -- Plotted selected ROI  --\n')

%% Overlay the subset with the WF image

pxl = 108;

cd('Z:\Christian-Sieben\data_HTP\2016-08-19_Nucleoid_MitoRNAGran\2016-08-19_Nuleoid_A647_WF4');
WF_name = '2016-08-19_Nuleoid_A647_WF4_MMStack_Pos0.ome.tif'
I = imread(WF_name);

CFX = (ceil((max(locs_Ch1(:,xCol)/pxl)))+1)./size(I);
CFY = (ceil((max(locs_Ch1(:,yCol)/pxl)))+1)./size(I);

figure('Position',[10 600 500 500],'name','Extracted Particles');
imshow(I,[1e2 2e3]); hold on;
scatter((subset1(minFrame:maxFrame,xCol)*CFX(:,1))/pxl,(subset1(minFrame:maxFrame,yCol)*CFY(:,1))/pxl,1,'green');hold on;
scatter((subset2(minFrame:maxFrame,xCol)*CFX(:,1))/pxl,(subset2(minFrame:maxFrame,yCol)*CFY(:,1))/pxl,1,'red');


%% Overlay all locs with the WF image

figure('Position',[10 600 500 500],'name','Extracted Particles');
imshow(I,[1e3 1e4]); hold on;
scatter(locs_Ch2(minFrame:maxFrame,xCol)/108,locs_Ch2(minFrame:maxFrame,yCol)/108,1,'red'); hold on;
scatter(locs_Ch1(minFrame:maxFrame,xCol)/108,locs_Ch1(minFrame:maxFrame,yCol)/108,1,'green'); hold on;

%% Create a control dataset for the respective FOV --> MALK Format

randVar = [];

randVar(:,1) = (xmax-xmin).*rand(length(subset1),1) + xmin;
randVar(:,2) = (ymax-ymin).*rand(length(subset1),1) + ymin;
randVar(:,3) = (max(subset1(:,frameCol))-min(subset1(:,frameCol))).*rand(length(subset1),1) + min(subset1(:,frameCol));
randVar(:,4) = (max(subset1(:,photonsCol))-max(subset1(:,photonsCol))).*rand(length(subset1),1) + min(subset1(:,photonsCol));

name_for_LAMA = [filename_locs_Ch1, '_ROI1_CBCcontrol_MALK.txt'];

path_Ch1 
fid = fopen(name_for_LAMA,'wt');
fprintf(fid, '# localization file (Malk format) \n# x[nm] \t	y[nm]	\t [frame] \t	I[a.u.]');

dlmwrite(name_for_LAMA, randVar,'delimiter', '\t', '-append')
fclose(fid);


fprintf('\n -- Control Dataset saved in Malk format --\n')


