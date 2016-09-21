clear, clc, close all
%% Load the raw data

filename_locs='2016-08-19_MitoRNAGran_A750_FOV_2_MMStack_Pos0_locResults_DC_corr'; 

name_with_header = '2016-08-19_MitoRNAGran_A750_FOV_2_MMStack_Pos0_locResults_DC';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

locs=dlmread([filename_locs '.dat'],',',1,0);
file = fopen([name_with_header '.dat']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol            = strmatch('x [nm]',h);
yCol            = strmatch('y [nm]',h);
dxCol           = strmatch('dx',h);
dyCol           = strmatch('dy',h);
frameCol        = strmatch('frame',h);
photonsCol      = strmatch('intensity [photon]',h);
sigmaCol        = strmatch('sigma [nm]',h);
loglikelihood   = strmatch('loglikelihood',h);
uncertainty     = strmatch('uncertainty [nm]',h); 

fprintf('\n -- Data loaded --\n')

%% Filter the data

% scatter(locs(:,xCol),locs(:,yCol),1)

histogram(locs(:,sigmaCol))

% Set Filter parameters 

minSigma            = 160;
maxSigma            = 130;
MinPhotons          = 200;
logFilter           = 75;
uncertainty_filter  = 20;
firstFrame          = 1000;

filter   = find(locs(:,frameCol) > firstFrame & locs(:,sigmaCol) < minSigma & locs(:,sigmaCol) > maxSigma & locs(:,photonsCol) > MinPhotons & locs(:,loglikelihood) < logFilter & locs(:,uncertainty) < uncertainty_filter );

subsetLL = [];
subsetLL = locs(filter,1:end);

scatter(subsetLL(:,xCol),subsetLL(:,yCol),1)

fprintf('\n -- Data Filtered (%f locs are left) --\n', ((length(subsetLL)/length(locs))));


%% Save Data in MALK format for LAMA


locs_subset (:,1) = subsetLL(:,xCol);
locs_subset (:,2) = subsetLL(:,yCol);
locs_subset (:,3) = subsetLL(:,frameCol);
locs_subset (:,4) = subsetLL(:,photonsCol);

name_for_LAMA = [filename_locs, '_MALK.txt'];

fid = fopen(name_for_LAMA,'wt');
fprintf(fid, '# localization file (Malk format) \n# x[nm] \t	y[nm]	\t [frame] \t	I[a.u.]');

dlmwrite(name_for_LAMA, locs_subset,'delimiter', '\t', '-append')
fclose(fid);


fprintf('\n -- Data saved in Malk format --\n')