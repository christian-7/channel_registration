clear, clc, close all


filename_locs_C1    ='2016-08-19_Nucleoid_A647_FOV_2_MMStack_Pos0_locResults_DC'; 

filename_locs_C2    ='2016-08-19_MitoRNAGran_A750_FOV_2_MMStack_Pos0_locResults_DC_corr'; 

name_with_header    ='2016-08-19_MitoRNAGran_A750_FOV_2_MMStack_Pos0_locResults_DC';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

locsC1=dlmread([filename_locs_C1 '.dat'],',',1,0);
locsC2=dlmread([filename_locs_C2 '.dat'],',',1,0);

file = fopen([name_with_header '.dat']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol = strmatch('x [nm]',h);
yCol = strmatch('y [nm]',h);
dxCol = strmatch('dx',h);
dyCol = strmatch('dy',h);
uncertainty = strmatch('uncertainty',h);
frameCol = strmatch('frame',h);
photonsCol = strmatch('intensity [photon]',h);


fprintf('\n -- Data loaded --\n')

%%  %%%%% Density filter %%%%%

NNdist=20;
idx = [];
idx = rangesearch(locsC1(:,xCol:yCol),locsC2(:,xCol:yCol),NNdist);

% search for locs in the A647 (C1) channel that are around locs in the A750 channel (C2)
% idx is a vector with length of locsC2, where for each loc a list of points in locsC1 is written
 
% clear all idx locs from A647 channel

locsC1_filtered = locsC1;
locsC1_filtered2 = [];

for i=1:length(idx);
    
locsC1_filtered(idx{i,1},1:end) = 0;    

end

for j = 1:12;
    
    if j==4; % Z column
    
        locsC1_filtered2(:,j) = 0;  
        
    else 
    locsC1_filtered2(:,j) = nonzeros(locsC1_filtered(:,j));
    
    end
    
end

length(locsC1_filtered2)/length(locsC1);
