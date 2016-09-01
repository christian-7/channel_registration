clear 
clc 
close all
%% Open data

filename_peaks='Au_fiducial_Roi2_2_PKF';

filename_peaks2=[filename_peaks '.prm'];
peaks=dlmread(filename_peaks2,',',1,0);

peaks(:,3)=peaks(:,3)*160;
peaks(:,4)=peaks(:,4)*160;

fprintf('\n -- Data loaded --\n')

%% Plot the data

figure
scatter(peaks(:,3),peaks(:,4),1,peaks(:,10));
colorbar

% scatter(peaks(:,20),peaks(:,21));

%% Select ROI
% 
xmin=min(peaks(:,3));
xmax=max(peaks(:,3));

ymin=min(peaks(:,4));
ymax=max(peaks(:,4));

% 
% xmin=3.82*1e4;
% xmax=3.86*1e4;
% 
% ymin=0.185*1e4;
% ymax=0.215*1e4;


vx=find(peaks(:,3)>xmin & peaks(:,3)<xmax);
subset1=peaks(vx,1:end);
vy=find(subset1(:,4)>ymin & subset1(:,4)<ymax);
subset2=subset1(vy,1:end);

figure
scatter(subset2(:,3),subset2(:,4))

%% Plot XY Scatter and XY drift

close all

% scatter(subset2(:,3)*100,subset2(:,4)*100,1);

scatter((subset2(:,3)-min(subset2(:,3))),(subset2(:,4)-min(subset2(:,4))));

figure('Position',[200 200 600 400],'Name','A549')
h=gcf;
set(h,'PaperOrientation','landscape');

subplot(2,1,1)
scatter(subset2(:,10),(subset2(:,3)-min(subset2(:,3))),1);
title('X drift')
xlabel('time [frame]')
ylabel('x position [nm]')

subplot(2,1,2)
scatter(subset2(:,10),(subset2(:,4)-min(subset2(:,4))),1);
title('Y drift')
xlabel('time [frame]')
ylabel('y position [nm]')

%% Fitting (figure for MS)

x=subset2(:,10);
y=subset2(:,3);

[fx,gofx,outputx] = fit(subset2(:,10),(subset2(:,3)-min(subset2(:,3))),'poly8'); % fit in nm, fx
[fy,gofy,outputy] = fit(subset2(:,10),(subset2(:,4)-min(subset2(:,4))),'poly8'); % fit in nm, fy

figure('Position',[200 200 400 600])
h=gcf;
set(h,'PaperOrientation','portrait');

subplot(2,1,1)
plot(fx,subset2(:,10),(subset2(:,3)-min(subset2(:,3))));
axis([0 max(peaks(:,10)) 0 700])
title('X drift')
xlabel('time [frames]')
ylabel('x position [nm]')
legend('data','fitted curve','Location','northeast');
legend('boxon');

subplot(2,1,2)
plot(fy,subset2(:,10),(subset2(:,4)-min(subset2(:,4))));
axis([0 max(peaks(:,10)) 0 700])
title('Y drift')
xlabel('time [frames]')
ylabel('y position [nm]')
legend('data','fitted curve','Location','northeast');
legend('boxon');

% Calculate RMSE

pdx=fitdist(outputx.residuals,'normal') % Distribution of residuals
pdy=fitdist(outputy.residuals,'normal') % Distribution of residuals


%% Calculate deviation for each frame

dx=[];
dy=[];

for frame=1:max(subset2(:,10));

dx(frame,1)=fx(frame)-fx(0); % in nm
dy(frame,1)=fy(frame)-fy(0); % in nm

end



%% Correct drift

subset2C=[];

j=1;

for i=1:length(subset2(:,10));
    
   
    
frame=subset2(i,10);

subset2C(j,1)=(subset2(i,3))+abs(dx(frame));
subset2C(j,2)=(subset2(i,4))-abs(dy(frame));

clear frame

j=j+1;

end
figure
scatter(subset2(:,3),subset2(:,4),2,'black');hold on;
scatter(subset2C(:,1),subset2C(:,2),2,'red');

%% Plot drift corrected bead

figure('Position',[200 200 400 600])
h=gcf;
set(h,'PaperOrientation','portrait');

subplot(2,1,1)
scatter(((subset2(:,3)-min(subset2(:,3))))+50,((subset2(:,4)-min(subset2(:,4))))+50,1,'black');
axis([0 400 0 400])
axis square 
box on
title('uncorrected')
xlabel('x position [nm]')
ylabel('y position [nm]')

subplot(2,1,2)
% scatter(((subset2C(:,1)-min(subset2C(:,1))))+50,((subset2C(:,2)-min(subset2C(:,2))))+50,1,'black');
scatter(((subset2C(:,1)-min(subset2(:,3))))+50,((subset2C(:,2)-min(subset2(:,4))))+50,1,'black');
% scatter(subset2C(:,1)+50,subset2C(:,2)+50,1,'black');
axis([0 400 0 400])
axis square 
box on
title('drift corrected')
xlabel('x position [nm]')
ylabel('y position [nm]')


pdCx=fitdist(subset2C(:,1),'normal') % Distribution of residuals
pdCy=fitdist(subset2C(:,2),'normal') % Distribution of residuals

%% Final figure

figure('Position',[200 200 600 600])
h=gcf;
set(h,'PaperOrientation','portrait');

subplot(2,2,1)
plot(fx,subset2(:,10),(subset2(:,3)-min(subset2(:,3))));
axis([0 max(peaks(:,10)) 0 500])
title('X drift')
xlabel('time [frames]')
ylabel('x position [nm]')
legend('data','fitted curve','Location','northeast');
legend('boxon');

subplot(2,2,3)
plot(fy,subset2(:,10),(subset2(:,4)-min(subset2(:,4))));
axis([0 max(peaks(:,10)) 0 500])
title('Y drift')
xlabel('time [frames]')
ylabel('y position [nm]')
legend('data','fitted curve','Location','northeast');
legend('boxon');

subplot(2,2,2)
scatter(((subset2(:,3)-min(subset2(:,3))))+50,((subset2(:,4)-min(subset2(:,4))))+50,1,'black');
axis([0 500 0 500])
axis square 
box on
title('uncorrected')
xlabel('x position [nm]')
ylabel('y position [nm]')

subplot(2,2,4)
% scatter(((subset2C(:,1)-min(subset2C(:,1))))+50,((subset2C(:,2)-min(subset2C(:,2))))+50,1,'black');
scatter(((subset2C(:,1)-min(subset2(:,3))))+50,((subset2C(:,2)-min(subset2(:,4))))+50,1,'black');
% scatter(subset2C(:,1)+50,subset2C(:,2)+50,1,'black');
axis([0 500 0 500])
axis square 
box on
title('drift corrected')
xlabel('x position [nm]')
ylabel('y position [nm]')

%% Open data

filename_peaks='Grouped5_with0gap_GPF';

filename_peaks2=[filename_peaks '.prm'];
peaks=dlmread(filename_peaks2,',',1,0);

peaks(:,3)=peaks(:,3)*160;
peaks(:,4)=peaks(:,4)*160;

fprintf('\n -- Data loaded --\n')

%% Calculate Sigma

% Sigma XY = 17 18
% Group Sigma XY = 22 23
% Group index = 19 
% Group size = 25 


figure('Position',[200 200 900 300])
h=gcf;
set(h,'PaperOrientation','landscape');

subplot(1,3,1)
bins=0:1:50;
xhist=hist(peaks(:,25),bins);
bar(bins,xhist/sum(xhist));
axis([0 20 0 0.5])
axis square 
box on
% title('Number of blinks')
xlabel('number of blinks')
ylabel('norm counts')

subplot(1,3,2)
bins=0:0.01:0.5;
shist=hist(peaks(:,17),bins);
bar(bins*160,shist/sum(shist));
axis([0 80 0 0.2])
axis square 
box on
title('before grouping')
xlabel('localization precision [nm]')
ylabel('norm counts')


subplot(1,3,3)
bins=0:0.01:0.5;
shist=hist(peaks(:,23),bins);
bar(bins*160,shist/sum(shist));
axis([0 80 0 0.2])
axis square 
box on
title('after grouping')
xlabel('localization precision [nm]')
ylabel('norm counts')


mean(peaks(:,23))*160
mean(peaks(:,17))*160


% hist(peaks(:,18),10);
% 
% pd_sigmax=fitdist(peaks(:,17),'normal')
% pd_sigmay=fitdist(peaks(:,18),'normal')


% hist(peaks(:,22),20);
% hist(peaks(:,23),20);
% 
% pd_sigmax=fitdist(peaks(:,22),'normal')
% pd_sigmay=fitdist(peaks(:,23),'normal')









