clear, clc, close all
%% Load the raw data

filename_locs = 'A549_SNA_A647_8_MMStack_Pos0_locResults_DC_TS'; 

name_with_header = 'A549_SNA_A647_8_MMStack_Pos0_locResults_cleaned';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

locs=dlmread([filename_locs '.dat'],',',1,0);
file = fopen([filename_locs '.dat']);
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


xCol            = strmatch('"x [nm]"',h);
yCol            = strmatch('"y [nm]"',h);
dxCol           = strmatch('dx',h);
dyCol           = strmatch('dy',h);
frameCol        = strmatch('"frame"',h);
photonsCol      = strmatch('"intensity [photon]"',h);
sigmaCol        = strmatch('"sigma [nm]"',h);
loglikelihood   = strmatch('"loglikelihood"',h);
uncertainty     = strmatch('"uncertainty [nm]"',h); 

fprintf('\n -- Data loaded --\n')

%% Filter the data

channel = 1;

% scatter(locs(:,xCol),locs(:,yCol),1)

histogram(locs(:,sigmaCol))

if channel==2;

% Set Filter parameters A750

minSigma            = 120;
maxSigma            = 250;
MinPhotons          = 200;
logFilter           = 150;
uncertainty_filter  = 20;
firstFrame          = 100;

else

% Set Filter parameters A647

maxSigma            = 260;
minSigma            = 80;
MinPhotons          = 500;
logFilter           = 75;
uncertainty_filter  = 20;
firstFrame          = 100;

end

filter   = find(locs(:,frameCol) > firstFrame & locs(:,sigmaCol) < maxSigma & locs(:,sigmaCol) > minSigma & locs(:,photonsCol) > MinPhotons & locs(:,loglikelihood) < logFilter & locs(:,uncertainty) < uncertainty_filter );

subsetLL = [];
subsetLL = locs(filter,1:end);

scatter(subsetLL(:,xCol),subsetLL(:,yCol),1)

fprintf('\n -- Data Filtered (%f locs are left) --\n', ((length(subsetLL)/length(locs))));


%% Save Data in MALK format for LAMA

locs_subset = [];

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

%% Create a control dataset for the respective FOV --> MALK Format

xmin = 3e4;
xmax = 5e4;
ymin = 3.4e4;
ymax = 5.1e4;

v1 = find(subsetLL(:,xCol) > xmin & subsetLL(:,xCol) < xmax & subsetLL(:,yCol) > ymin & subsetLL(:,yCol) < ymax);

subset1=subsetLL(v1,1:end);

figure('Position',[100 400 500 500],'Name','Drift Ch1')
scatter(subset1(:,xCol),subset1(:,yCol),1,'bo')

fprintf('\n -- Plotted selected ROI  --\n')

randVar = [];

randVar(:,1) = (xmax-xmin).*rand(length(subset1),1) + xmin;
randVar(:,2) = (ymax-ymin).*rand(length(subset1),1) + ymin;
randVar(:,3) = (max(subset1(:,frameCol))-min(subset1(:,frameCol))).*rand(length(subset1),1) + min(subset1(:,frameCol));
randVar(:,4) = (max(subset1(:,photonsCol))-max(subset1(:,photonsCol))).*rand(length(subset1),1) + min(subset1(:,photonsCol));

name_for_LAMA = [filename_locs, '_CBC_control_MALK.txt'];
 
fid = fopen(name_for_LAMA,'wt');
fprintf(fid, '# localization file (Malk format) \n# x[nm] \t	y[nm]	\t [frame] \t	I[a.u.]');

dlmwrite(name_for_LAMA, randVar,'delimiter', '\t', '-append')
fclose(fid);

fprintf('\n -- Control Dataset saved in Malk format --\n')

