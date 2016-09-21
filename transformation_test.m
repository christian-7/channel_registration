%  Load the bead coordinates from both channels

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

%% 

tform = estimateGeometricTransform(moving,fixed,'affine')

[x,y] = transformPointsForward(tform,moving(:,1), moving(:,2))

figure;hold all;
scatter(fixed(:,1),fixed(:,2),'gx');
scatter(moving(:,1),moving(:,2),'ro')
scatter(x,y,'rx')

%% second order polynomial warp transform
% 
% fixed = [1.5 4; 3 2; 5 3; 4 4; 5 5; 6 6];
% moving = [2 3; 3.5 4; 4 5; 5 6; 6 7; 7 8];
% 

% First transformation

Tgr = cp2tform(fixed, moving, 'similarity');
[xgt, ygt] = tformfwd(Tgr, fixed(:,1), fixed(:,2));
afterFirstTF=[xgt ygt];

% Second transformation
% must be applied to the NN in the dataset

% T_lwm = cp2tform(afterFirstTF,moving,'lwm',36);'polynomial',4

T_lwm = cp2tform(afterFirstTF,moving,'polynomial',2);
[xgt2, ygt2] = tformfwd(T_lwm, xgt, ygt);

% tform = images.geotrans.PolynomialTransformation2D(moving,fixed,2)
 
% movingPointsEstimated = transformPointsInverse(tform,fixed);

figure
scatter(fixed(:,1),fixed(:,2),5,'green');hold on;
scatter(moving(:,1),moving(:,2),5,'red');hold on;
scatter(xgt,ygt,5,'+');hold on;
scatter(movingPointsEstimated(:,1),movingPointsEstimated(:,2),'*');hold on;


errorInFit = hypot(movingPointsEstimated(:,1)-moving(:,1),...
                   movingPointsEstimated(:,2)-moving(:,2));
               
               
               

%% Rigid transform



% find same bead in both channels
% calculate dX and dY for each frame
% apply to dataset



