%% Calculate 2D affine tranformation for 2 sets of fiducial coordinates

% Input : fiducial coordinates in TS format
% Output: affine_transform.mat to use with apply_affine_transform.m

% Load the bead coordinates from both channels

close all, clear, clc

% The script calculates the correction for the A750 channel

cd('.\test_data');

filename_peaksc1 = 'Beads_Ch_A647';             % filename of TS output file
filename_peaksc2 = 'Beads_Ch_A750';             % filename of TS output file

% filename_peaksc1='Ch1_locResult';             % filename of TS output file
% filename_peaksc2='Ch2_locResult';             % filename of TS output file

peaksC1=dlmread([filename_peaksc1 '.csv'],',',1,0);
peaksC2=dlmread([filename_peaksc2 '.csv'],',',1,0);

file = fopen([filename_peaksc1 '.csv']);
line = fgetl(file);
h    = regexp( line, ',', 'split' );

xCol = strmatch('"x [nm]"',h);
yCol = strmatch('"y [nm]"',h);

fixed  = peaksC1(:,2:3);
moving = peaksC2(:,2:3); % The script calculates the correction for the A750 channel

cd('.\..')

% Filter only points that appear in both channels

Idx_c1 = rangesearch(moving,fixed,500);
Idx_c2 = rangesearch(fixed,moving,500);

selection_c1 = []; selection_c2 = [];

for i = 1:length(Idx_c1);
    
    if isempty(Idx_c1{i});
        
    else
    
    selection_c1 = vertcat(selection_c1, i);    
    
    end
    
end


for i = 1:length(Idx_c2);
    
    if isempty(Idx_c2{i});
        
    else
    
    selection_c2 = vertcat(selection_c2, i);    
    
    end
    
end


fixed_s = []; moving_s = [];

fixed_s   = fixed(selection_c1,1:end);
moving_s  = moving(selection_c2,1:end);

figure('Position',[800 500 600 500]); hold all;
scatter(fixed_s(:,1),fixed_s(:,2),'rx');
scatter(moving_s(:,1),moving_s(:,2),'gx');  
legend('seleclted points Ch1','seleclted points Ch2');
title('2C data before transformation');
box on;
xlabel('x(nm)')
ylabel('y(nm)')


% Calculate and apply transformation

tform = estimateGeometricTransform(moving_s,fixed_s,'affine');
[corrected(:,1),corrected(:,2)] = transformPointsForward(tform,moving_s(:,1), moving_s(:,2));

% tform       = fitgeotrans(moving_s,fixed_s,'polynomial',2);
% corrected   = transformPointsInverse(tform,fixed_s);

figure('Position',[100 500 600 500]); hold all;
scatter(fixed_s(:,1),fixed_s(:,2),'gx');
scatter(moving_s(:,1),moving_s(:,2),'ro')
scatter(corrected(:,1),corrected(:,2),'rx')
box on;
title('2C data after transformation');
legend('Channel A750','Channel A647','Channel A647 corrected')
xlabel('x(nm)')
ylabel('y(nm)')

deviation_error = [];

Idx_error = rangesearch(fixed_s,corrected,200);

for i = 1:length(Idx_error);
    
    if isempty(Idx_error{i});
        
    else
      
    deviation_error(i,1) = sqrt(((corrected(i,1) - fixed_s(Idx_error{i},1))^2) + ((corrected(i,2) - fixed_s(Idx_error{i},2))^2));
    
    end
end

meanError = mean(deviation_error)

stdError = std(deviation_error)

figure('Position',[800 100 600 500])
scatter(corrected(:,1),corrected(:,2),15,deviation_error,'filled')
colorbar
title(['Color -> registration error, Mean = ', num2str(meanError)]);
box on
xlabel('x(nm)')
ylabel('y(nm)')


save('aff_transform','tform');

fprintf('\n -- Transformation saved --\n')


