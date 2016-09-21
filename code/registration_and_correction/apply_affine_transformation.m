%% Applies both the affine and the linear transform

% Input:        A750 localization file, drift-corrected
%               dx, dy
%               affine_transform

% Output:       corrected A750 localization file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the bead coordinates from both channels

close all, clear, clc

load('D:\Christian\GitHub\channel_registration\test_data\aff_transform.mat')

% The script calculates the correction for the A750 channel

filename_peaksc1 = '2016-08-19_MitoRNAGran_A750_FOV_5_MMStack_Pos0_locResults_DC';             % filename of TS output file

peaks=dlmread([filename_peaksc1 '.dat'],',',1,0);

file = fopen([filename_peaksc1 '.dat']);
line = fgetl(file);
h    = regexp( line, ',', 'split' );

xCol        = strmatch('x [nm]',h);
yCol        = strmatch('y [nm]',h);
frameCol    = strmatch('frame',h);

fprintf('\n -- Data loaded --\n') 

%% Apply the transformation

[corrected(:,1),corrected(:,2)] = transformPointsForward(tform,peaks(:,xCol), peaks(:,yCol));

peaks2 = peaks;
peaks2(:,xCol) = corrected(:,1);
peaks2(:,yCol) = corrected(:,2);

% filename_peaksc2 = [filename_peaksc1 '_transformed.dat'];
% dlmwrite(filename_peaksc2,peaks2);

fprintf('\n -- Transformation applied --\n') 

%% Apply linear transformation

corr_Ch2 = peaks2;

load('devX.mat')
load('devY.mat')

for i=1:length(corr_Ch2(:,frameCol));
    
frame = corr_Ch2(i,frameCol);

corr_Ch2(i,xCol) = (corr_Ch2(i,xCol))-abs(dx(frame)-50);
corr_Ch2(i,yCol) = (corr_Ch2(i,yCol))+abs(dy(frame)+50);

clear frame

end

dlmwrite([filename_peaksc1 '_corr.dat'],corr_Ch2);

fprintf('\n -- Saved corrected localization file --\n') 


