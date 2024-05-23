clear;
clc;

inpath = './seismograms/';

[seisx, nt, nx, dt] = readsu([inpath, 'seisxFD_S.su']);
[seisz] = readsu([inpath, 'seiszFD_S.su']);


% define colormap
N = 1000;
cmap = ([linspace(180, 255, N)', linspace(28 , 255, N)', linspace(27 , 255, N)'; ...
         linspace(254,   5, N)', linspace(254,   4, N)', linspace(254, 131, N)'])/255;


x = (1:1:nx)';
t = (0:1:nt-1)'*dt;


perc = 1.0;

figure(1);
subplot(2,1,1);
zmin = min(seisx(:));
zmax = max(seisx(:));
zclip = min(-zmin, zmax);
hold off;
imagesc(x, t, seisx, perc*[-zclip, zclip]);
hold off;
colormap(gca, cmap);
xlabel('Trace#');
ylabel('Time [s]');
title('Ridial');

subplot(2,1,2);
zmin = min(seisz(:));
zmax = max(seisz(:));
zclip = min(-zmin, zmax);
hold off;
imagesc(x, t, seisz, perc*[-zclip, zclip]);
hold off;
colormap(gca, cmap);
xlabel('Trace#');
ylabel('Time [s]');
title('Vertical');
