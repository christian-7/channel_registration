%% Two color channel registration using local weighted mean transformation 

% Reference: L.S. Churchman, J.A. Spudich, Single-molecule high-resolution colocalization of single probes, Cold Spring Harb. Protoc. 2012 (2012) 242–245

% Input: 

%           bead coordinates in both channels
%           localization files 

% Output: 

%           transformation function
%           corrected dataset


%% Load bead data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option1: Load the bead data % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('.\test_data');

filename_beads_642='Ch1_locResult';             % filename of TS output file
filename_beads_750='Ch2_locResult';             % filename of TS output file

peaksC1=dlmread([filename_beads_642 '.csv'],',',1,0);
peaksC2=dlmread([filename_beads_750 '.csv'],',',1,0);

file = fopen([filename_beads_642 '.csv']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol = strmatch('"x [nm]"',h);
yCol = strmatch('"y [nm]"',h);

fixed = peaksC1(:,2:3);
moving = peaksC2(:,2:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option2:  Load the localization data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load from HDF input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option3:  Load the WF images of beads taken in each Channel %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('.\test_data');

Ch1 = imread('A647_1.tif');
Ch2 = imread('A750_1.tif');

rgbIm = cat(3,Ch1,Ch2,zeros(size(Ch1)));
imshow(rgbIm);

[Ch1_Pts,Ch2_Pts]= cpselect(Ch1, Ch2,'Wait',true);

fixed = Ch1_Pts;
moving_1 = Ch2_Pts;

%% Calculate LWM transformation

[IDX] = rangesearch(moving,fixed,500); % find the point in moving that is closest to the point in fixed

for i = 1:length(IDX);

    moving_1(i,1) = moving(IDX{i,1},1);
    moving_1(i,2) = moving(IDX{i,1},2);

end

T_lwm = fitgeotrans(fixed,moving_1,'lwm',10);

corrected_moving = transformPointsInverse(T_lwm,moving_1);


for i = 1:length(fixed);
    
    TRE(:,i) = sqrt((fixed(i,1)-corrected_moving(i,1))^2 + (fixed(i,1)-corrected_moving(i,1))^2);

end

figure('Position',[100 600 500 500])

scatter(fixed(:,1),fixed(:,2),'gx');hold on;
scatter(moving_1(:,1),moving_1(:,2),'ro');hold on;
scatter(corrected_moving(:,1),corrected_moving(:,2),'rx');hold on;
legend('Fixed', 'Moving','Corrected');
title(['TRE = ' num2str(mean(TRE))]);
box on; axis square


diff = [];
diff = fixed - moving_1;


figure('Position',[700 600 500 500])

scatter(fixed(:,1),fixed(:,2),'gx');hold on;
scatter(moving_1(:,1),moving_1(:,2),'ro');hold on;

for i = 1:length(fixed);

q = quiver(moving_1(i,1),moving_1(i,2),diff(i,1),diff(i,2)); hold on;
q.Color = 'black'
q.LineWidth = 2;
q.MaxHeadSize = 2;

box on; axis square

title(['TRE = ' num2str(mean(TRE))]);
end

%% Apply transformation to A750 dataset
