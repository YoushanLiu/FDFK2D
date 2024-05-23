clear;
clc;

inpath = './input/';

[vp, nz, nx] = readsu([inpath, 'vp_model_smooth_2.4x1.2_without_pml.su']);
vs = readsu([inpath, 'vs_model_smooth_2.4x1.2_without_pml.su']);


x0 = 0.0;
z0 = 0.0;
dx = 200.0;
dz = 200.0;

x = (x0 + (nx-1)*dx)*1.e-3;
z = (z0 + (nz-1)*dz)*1.e-3;


figure(1);
subplot(2,1,1);
hold off;
imagesc(x, z, vp*1.e-3);
hold off;
xlabel('Distance [km]');
ylabel('Depth [km]');
title('Vp');


subplot(2,1,2);
hold off;
imagesc(x, z, vs*1.e-3);
hold off;
xlabel('Distance [km]');
ylabel('Depth [km]');
title('Vs');