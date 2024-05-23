clear;
clc;

inpath = './seismograms/';

%[ux, nt, nx, dt] = readsu([inpath, 'seisx.su']);
%[uz] = readsu([inpath, 'seisz.su']);

inpath = './seismograms_1layer_P_S/';

[ux, nt, nx, dt] = readsu([inpath, 'seisxFD_P.su']);
[uz] = readsu([inpath, 'seiszFD_P.su']);


itr = 101;

Tmax = (nt-1)*dt;

t = (0:1:nt-1)'*dt;


% parameters for RFs
tlag_min = -2.0;
tlag_max = 0.5*Tmax + tlag_min;
tlag = (tlag_min:dt:tlag_max);


% compute RF
RF = makeRFitdecon(ux(:,itr), uz(:,itr), dt, tlag_min, tlag_max, 0.0, 8., 100, 1.e-3, 1);


figure(1);
subplot(3,1,1);
hold off;
plot(t, ux(:,itr), 'r');
hold off;
xlabel('Time [s]');
ylabel('Amplitude');
title('Radial');

subplot(3,1,2);
hold off;
plot(t, uz(:,itr), 'r');
hold off;
xlabel('Time [s]');
ylabel('Amplitude');
title('Vertical');

subplot(3,1,3);
hold off;
plot(tlag, RF, 'r');
hold off;
xlabel('Lag [s]');
ylabel('Amplitude');
title('RF');