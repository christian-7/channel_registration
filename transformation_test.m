%% Test different point transformations

% affine, polynomial, lwm

% LWM from : L.S. Churchman, J.A. Spudich, Single-molecule high-resolution colocalization of single probes, Cold Spring Harb. Protoc. 2012 (2012) 242�245

clear, clc, close all

%%  Load the test data. Bead coordinates from both channels

cd('.\test_data');

filename_peaksc1='Ch1_locResult';             % filename of TS output file
filename_peaksc2='Ch2_locResult';             % filename of TS output file

peaksC1=dlmread([filename_peaksc1 '.csv'],',',1,0);
peaksC2=dlmread([filename_peaksc2 '.csv'],',',1,0);

file = fopen([filename_peaksc1 '.csv']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol = strmatch('"x [nm]"',h);
yCol = strmatch('"y [nm]"',h);

fixed = peaksC1(:,2:3);
moving = peaksC2(:,2:3);

%% Select Points from reference image

cd('.\test_data');

Ch1 = imread('A647_1.tif');
Ch2 = imread('A750_1.tif');

rgbIm = cat(3,Ch1,Ch2,zeros(size(Ch1)));
imshow(rgbIm);

[Ch1_Pts,Ch2_Pts]= cpselect(Ch1, Ch2,'Wait',true);

fixed = Ch1_Pts;
moving_1 = Ch2_Pts;

%% Affine transformation

[IDX] = rangesearch(moving,fixed,500); % find the point in moving that is closest to the point in fixed

for i = 1:length(IDX);

    moving_1(i,1) = moving(IDX{i,1},1);
    moving_1(i,2) = moving(IDX{i,1},2);

end

tform = estimateGeometricTransform(moving_1,fixed,'affine');

[corrected_moving(:,1),corrected_moving(:,2)] = transformPointsForward(tform,moving_1(:,1), moving_1(:,2));


for i = 1:length(fixed);
    
    TRE(:,i) = sqrt((fixed(i,1)-corrected_moving(i,1))^2 + (fixed(i,1)-corrected_moving(i,1))^2);

end

figure
scatter(fixed(:,1),fixed(:,2),'gx');hold on;
scatter(moving_1(:,1),moving_1(:,2),'ro');hold on;
scatter(corrected_moving(:,1),corrected_moving(:,2),'rx');hold on;
legend('Fixed', 'Moving','Corrected');
title(['TRE = ' num2str(mean(TRE))]);

%% Using a polynomial transformation

[IDX] = rangesearch(moving,fixed,500); % find the point in moving that is closest to the point in fixed

for i = 1:length(IDX);

    moving_1(i,1) = moving(IDX{i,1},1);
    moving_1(i,2) = moving(IDX{i,1},2);

end

T_poly = fitgeotrans(fixed,moving_1,'polynomial',2);
corrected_moving = transformPointsInverse(T_poly,moving_1);

for i = 1:length(fixed);
    
    TRE(:,i) = sqrt((fixed(i,1)-corrected_moving(i,1))^2 + (fixed(i,1)-corrected_moving(i,1))^2);

end

figure
scatter(fixed(:,1),fixed(:,2),'gx');hold on;
scatter(moving_1(:,1),moving_1(:,2),'ro');hold on;
scatter(corrected_moving(:,1),corrected_moving(:,2),'rx');hold on;
legend('Fixed', 'Moving','Corrected');
title(['TRE = ' num2str(mean(TRE))]);

%% Using a local weihgted mean transformation

[IDX] = rangesearch(moving,fixed,500); % find the point in moving that is closest to the point in fixed

for i = 1:length(IDX);

    moving_1(i,1) = moving(IDX{i,1},1);
    moving_1(i,2) = moving(IDX{i,1},2);

end

T_lwm = fitgeotrans(fixed,moving_1,'lwm',7);

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
q.LineWidth = 1;
q.MaxHeadSize = 2;

box on; axis square

title(['TRE = ' num2str(mean(TRE))]);
end

