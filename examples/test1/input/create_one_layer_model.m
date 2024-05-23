clear;
clc;

npml = 20;
X = 200.e3;
Z = 50.e3;

dx = 200;
dz = 200;


nx = round(X/dx)+1 + 2*npml;
nz = round(Z/dz)+1 + npml;

rho = zeros(nz,nx);
vp  = zeros(nz,nx);
vs  = zeros(nz,nx);

rho1 = 2700.0; Vp1 = 6000.0; Vs1 = 3450.0;
rho2 = 3300.0; Vp2 = 8000.0; Vs2 = 4480.0;


h = 35.e3;
ih = round(h/dz);


rho(:,:) = rho1;
vp(:,:)  = Vp1;
vs(:,:)  = Vs1;


figure(1);
subplot(3,1,1);
hold off;
imagesc(rho);
hold off;
xlabel('Distance [km]');
ylabel('Depth [km]');
title('rho');

subplot(3,1,2);
hold off;
imagesc(vp);
hold off;
xlabel('Distance [km]');
ylabel('Depth [km]');
title('Vp');

subplot(3,1,3);
hold off;
imagesc(vs);
hold off;
xlabel('Distance [km]');
ylabel('Depth [km]');
title('Vs');


writesu('rho_model.su', rho, dz, 1);
writesu('vp_model.su', vp, dz, 1);
writesu('vs_model.su', vs, dz, 1);