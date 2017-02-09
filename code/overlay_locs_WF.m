%% Use before CBC analysis 

% What this does:

% Visualize the FOV 
% choose ROI for CBC analysis
% export overlay images
% create a control dataset

clear, clc, close all
%% Load the two datasets 

path_Ch1            = 'Z:\Christian-Sieben\data_HTP\2016-09-22_A549_EGFR_SNA\locResults\A549_SNA_A647_11';
filename_locs_Ch1   = 'A549_SNA_A647_11_MMStack_Pos0_locResults_DC_TS'; 

path_Ch2            = 'Z:\Christian-Sieben\data_HTP\2016-09-22_A549_EGFR_SNA\locResults\A549_EGFR_A750_11';
filename_locs_Ch2   = 'A549_EGFR_A750_11_MMStack_Pos0_locResults_DC_corr'; 

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
frameCol    = strmatch('frame',h);
photonsCol  = strmatch('intensity [photon]',h);

xCol        = strmatch('"x [nm]"',h);
yCol        = strmatch('"y [nm]"',h);
frameCol    = strmatch('"frame"',h);
photonsCol  = strmatch('"intensity [photon]"',h);


fprintf('\n -- Data loaded --\n')

%% Show 2 color dataset

minFrame = 1;
maxFrame = 1e4;

figure; hold all;
scatter(locs_Ch1(minFrame:maxFrame,xCol),locs_Ch1(minFrame:maxFrame,yCol),1,'green');
scatter(locs_Ch2(minFrame:maxFrame,xCol),locs_Ch2(minFrame:maxFrame,yCol),1,'red');

%% Select ROI

xmin = 3.6e4;
xmax = 5.2e4;
ymin = 3.4e4;
ymax = 4.8e4;


v1 = find(locs_Ch1(:,xCol) > xmin & locs_Ch1(:,xCol) < xmax & locs_Ch1(:,yCol) > ymin & locs_Ch1(:,yCol) < ymax);
v2 = find(locs_Ch2(:,xCol) > xmin & locs_Ch2(:,xCol) < xmax & locs_Ch2(:,yCol) > ymin & locs_Ch2(:,yCol) < ymax);

subset1=locs_Ch1(v1,1:end);
subset2=locs_Ch2(v2,1:end);

figure('Position',[100 400 500 500],'Name','Drift Ch1')
scatter(subset1(minFrame:maxFrame,xCol),subset1(minFrame:maxFrame,yCol),1,'bo')

figure('Position',[1200 400 500 500],'Name','Drift Ch2')
scatter(subset2(minFrame:maxFrame,xCol),subset2(minFrame:maxFrame,yCol),1,'bo')

figure('Position',[1200 400 500 500],'Name','Overlay')
scatter(subset1(minFrame:maxFrame,xCol),subset1(minFrame:maxFrame,yCol),1,'green'); hold on;
scatter(subset2(minFrame:maxFrame,xCol),subset2(minFrame:maxFrame,yCol),1,'red')

fprintf('\n -- Plotted selected ROI  --\n')

%% Overlay the subset with the WF image

pxl = 108;

cd('Z:\Christian-Sieben\data_HTP\2016-08-19_Nucleoid_MitoRNAGran\2016-08-19_Nuleoid_A647_WF2');
WF_name = '2016-08-19_Nuleoid_A647_WF2_MMStack_Pos0.ome.tif';
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


%% Render the ROI and save with Ch1 data

pxlsize = 10;

size(1,1)    = round((max(subset1(:,yCol))-min(subset1(:,yCol)))/pxlsize);
size(1,2)    = round((max(subset1(:,xCol))-min(subset1(:,xCol)))/pxlsize);

size(2,1)    = round((max(subset2(:,yCol))-min(subset2(:,yCol)))/pxlsize);
size(2,2)    = round((max(subset2(:,xCol))-min(subset2(:,xCol)))/pxlsize);

im1 = hist3([subset1(:,yCol),subset1(:,xCol)],[min(size(:,1)) min(size(:,2))]);
im2 = hist3([subset2(:,yCol),subset2(:,xCol)],[min(size(:,1)) min(size(:,2))]);

imG1 = uint16(imgaussfilt(im1, pxlsize/10));
imG2 = uint16(imgaussfilt(im2, pxlsize/10));

cd(path_Ch1);


imwrite(im1,'imageCh1_A647.tiff');
imwrite(im2,'imageCh2_A750.tiff');

imwrite(imG1,'imageCh1_A647_16bit.tiff');
imwrite(imG2,'imageCh2_A750_16bit.tiff');
% 
% name = [filename_locs_Ch1 '_rendered_' num2str(pxlsize) 'nm_per_pxl.tiff'];  
% 
% I32 = [];
% I32 = uint32(imG1_16b);
% 
% t = Tiff('name.tiff','w');
% 
% tagstruct.ImageLength     = size(I32,1);
% tagstruct.ImageWidth      = size(I32,2);
% tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
% tagstruct.BitsPerSample   = 32;
% tagstruct.SamplesPerPixel = 1;
% tagstruct.RowsPerStrip    = 16;
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% tagstruct.Software        = 'MATLAB';
% t.setTag(tagstruct)
% 
% t.write(I32);
% t.close()
% 
% 
% name = [filename_locs_Ch2 '_rendered_' num2str(pxlsize) 'nm_per_pxl.tiff'];  
% name = [savename '_rendered_' num2str(pxlsize) 'nm_per_pxl.tiff'];  
% 
% I32=[];
% I32=uint32(imG2);
% 
% t = Tiff(name,'w');
% 
% tagstruct.ImageLength     = size(I32,1);
% tagstruct.ImageWidth      = size(I32,2);
% tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
% tagstruct.BitsPerSample   = 32;
% tagstruct.SamplesPerPixel = 1;
% tagstruct.RowsPerStrip    = 16;
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% tagstruct.Software        = 'MATLAB';
% t.setTag(tagstruct)
% 
% t.write(I32);
% t.close()


% cd(locpath);

fprintf('\n -- Image saved --\n');


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


