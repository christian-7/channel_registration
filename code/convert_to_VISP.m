clear, clc, close all

filename_locs       = '2016-08-19_Nucleoid_A647_FOV_6_MMStack_Pos0_locResults_DC'; 

name_with_header    = '2016-08-19_Nucleoid_A647_FOV_6_MMStack_Pos0_locResults_DC';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

locs=dlmread([filename_locs '.dat'],',',1,0);
file = fopen([name_with_header '.dat']);
line = fgetl(file);
h = regexp( line, ',', 'split' );

xCol            = strmatch('x [nm]',h);
yCol            = strmatch('y [nm]',h);
dxCol           = strmatch('dx',h);
dyCol           = strmatch('dy',h);
uncertainty     = strmatch('uncertainty',h);
frameCol        = strmatch('frame',h);
photonsCol      = strmatch('intensity [photon]',h);
sigmaCol        = strmatch('sigma [nm]',h);
loglikelihood   = strmatch('loglikelihood',h);

% xCol            = strmatch('"x [nm]"',h);
% yCol            = strmatch('"y [nm]"',h);
% dxCol           = strmatch('dx',h);
% dyCol           = strmatch('dy',h);
% frameCol        = strmatch('"frame"',h);
% photonsCol      = strmatch('"intensity [photon]"',h);
% sigmaCol        = strmatch('"sigma [nm]"',h);
% loglikelihood   = strmatch('"loglikelihood"',h);
% uncertainty     = strmatch('"uncertainty [nm]"',h); 



fprintf('\n -- Data loaded --\n')

%% Filter the data

channel = 1;

% scatter(locs(:,xCol),locs(:,yCol),1)

histogram(locs(:,uncertainty))

if channel==2;

% Set Filter parameters A750

minSigma            = 160;
maxSigma            = 260;
MinPhotons          = 200;
logFilter           = 150;
uncertainty_filter  = 20;
firstFrame          = 100;

else

% Set Filter parameters A647

maxSigma            = 150;
minSigma            = 120;
MinPhotons          = 1500;
logFilter           = 75;
uncertainty_filter  = 30;
firstFrame          = 500;

end

filter   = find(locs(:,frameCol) > firstFrame & locs(:,sigmaCol) < maxSigma & locs(:,sigmaCol) > minSigma & locs(:,photonsCol) > MinPhotons & locs(:,loglikelihood) < logFilter & locs(:,uncertainty) < uncertainty_filter );

subsetLL = [];
subsetLL = locs(filter,1:end);

scatter(subsetLL(:,xCol),subsetLL(:,yCol),1)

fprintf('\n -- Data Filtered (%f locs are left) --\n', ((length(subsetLL)/length(locs))));

%% 

% VISP format x,y,cPrec, yPrec, Photons, Frame

locs_subset = [];

locs_subset (:,1) = subsetLL(:,xCol);
locs_subset (:,2) = subsetLL(:,yCol);
locs_subset (:,3) = subsetLL(:,uncertainty);
locs_subset (:,4) = subsetLL(:,uncertainty);
locs_subset (:,5) = subsetLL(:,photonsCol);
locs_subset (:,6) = subsetLL(:,frameCol);

name_for_VISP = [filename_locs, '_forVISP.txt'];

% fid = fopen(name_for_LAMA,'wt');
% fprintf(fid, '# localization file (Malk format) \n# x[nm] \t	y[nm]	\t [frame] \t	I[a.u.]');

dlmwrite(name_for_VISP, locs_subset,'delimiter', '\t');
% fclose(fid);


fprintf('\n -- Data saved in VISP format --\n')