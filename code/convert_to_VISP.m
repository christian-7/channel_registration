clear, clc, close all


filename_locs='2016-08-19_MitoRNAGran_A750_FOV_2_MMStack_Pos0_locResults_DC_corr'; 

name_with_header = '2016-08-19_MitoRNAGran_A750_FOV_2_MMStack_Pos0_locResults_DC';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

locs=dlmread([filename_locs '.dat'],',',1,0);
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


% VISP format x,y,cPrec, yPrec, Photons, Frame

locs_subset (:,1) = locs(:,xCol);
locs_subset (:,2) = locs(:,yCol);
locs_subset (:,3) = locs(:,uncertainty);
locs_subset (:,4) = locs(:,uncertainty);
locs_subset (:,5) = locs(:,photonsCol);
locs_subset (:,6) = locs(:,frameCol);

name_for_VISP = [filename_locs, '_forVISP.txt'];

% fid = fopen(name_for_LAMA,'wt');
% fprintf(fid, '# localization file (Malk format) \n# x[nm] \t	y[nm]	\t [frame] \t	I[a.u.]');

dlmwrite(name_for_VISP, locs_subset,'delimiter', '\t');
% fclose(fid);


fprintf('\n -- Data saved in VISP format --\n')