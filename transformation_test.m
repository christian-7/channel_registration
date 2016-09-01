%% second order polynomial warp transform

fixed = [1.5 4; 3 2; 5 3; 4 4; 5 5; 6 6];
moving = [2 3; 3.5 4; 4 5; 5 6; 6 7; 7 8];

tform = images.geotrans.PolynomialTransformation2D(moving,fixed,2)
 
movingPointsEstimated = transformPointsInverse(tform,fixed);

figure
scatter(fixed(:,1),fixed(:,2));hold on;
scatter(moving(:,1),moving(:,2),5,'+');hold on;
scatter(movingPointsEstimated(:,1),movingPointsEstimated(:,2),'*');hold on;

errorInFit = hypot(movingPointsEstimated(:,1)-moving(:,1),...
                   movingPointsEstimated(:,2)-moving(:,2))
               
%% Rigid transform



% find same bead in both channels
% calculate dX and dY for each frame
% apply to dataset



