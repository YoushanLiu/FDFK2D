clear;
clc;

isnap = 1000;

inpath = './snapshots/';

[ux, nz, nx] = readsu([inpath, sprintf('%5.5d', isnap), 'ux.su']);
[uz] = readsu([inpath, sprintf('%5.5d', isnap), 'uz.su']);


% define colormap
N = 1000;
cmap = ([linspace(180, 255, N)', linspace(28 , 255, N)', linspace(27 , 255, N)'; ...
         linspace(254,   5, N)', linspace(254,   4, N)', linspace(254, 131, N)'])/255;



x0 = 0.0;
z0 = 0.0;
dx = 200.0;
dz = 200.0;

x = (x0 + (nx-1)*dx)*1.e-3;
z = (z0 + (nz-1)*dz)*1.e-3;


perc = 1.0;

figure(1);
subplot(2,1,1);
zmin = min(ux(:));
zmax = max(ux(:));
zclip = min(-zmin, zmax);
hold off;
imagesc(x, z, ux, perc*[-zclip, zclip]);
hold off;
colormap(gca, cmap);
xlabel('Distance [km]');
ylabel('Depth [km]');
title('Horizontal');

subplot(2,1,2);
zmin = min(uz(:));
zmax = max(uz(:));
zclip = min(-zmin, zmax);
hold off;
imagesc(x, z, uz, perc*[-zclip, zclip]);
hold off;
colormap(gca, cmap);
xlabel('Distance [km]');
ylabel('Depth [km]');
title('Vertical');

