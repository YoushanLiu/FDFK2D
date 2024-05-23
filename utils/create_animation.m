clear;
clc;


% define colormap
N = 1000;
cmap = ([linspace(180, 255, N)', linspace(28 , 255, N)', linspace(27 , 255, N)'; ...
         linspace(254,   5, N)', linspace(254,   4, N)', linspace(254, 131, N)'])/255;



nsnap = 800;

nx = 2401;
nz = 601;

dx = 0.2;
dz = 0.2;

zMoho = 54;

dt = 0.01;

xo = -50;
zo = 0.0;
x = xo + (0:1:nx-1)*dx;
z = zo + (0:1:nz-1)*dz;

x0 = min(x); x1 = max(x);
z0 = min(z); z1 = max(z);
X = (x1-x0);
Z = (z1-z0);


% % abnormal bodies
% Blength = 20;
% xb0 = x0 + X/2 - Blength/2;
% xb1 = x0 + X/2 + Blength/2;
% 
% % Bheight = 8;
% % Bdepth = 8;
% % zb0 = z0 + Bdepth - Bheight/2;
% % zb1 = z0 + Bdepth + Bheight/2;
% Btop = 5;
% Bheight = 8;
% zb0 = z0 + Btop;
% zb1 = z0 + Btop + Bheight;
% 
% 
% Blength = 150;
% x_reserved = 15;
% xo_left  = x_reserved;
% xo_right = Blength - x_reserved;
% x_transition_zone = 25;
% xb3 = xo_left  + x_transition_zone;
% xb4 = xo_right - x_transition_zone;
% Bdepth_bottom = 2;
% zb3 = 0;
% zb4 = zb3 + Bdepth_bottom;


% % read boundaries
% line1 = load('Line1.txt', '-ascii');
% line2 = load('Line2.txt', '-ascii');
% line3 = load('Line3.txt', '-ascii');
% line4 = load('Line4.txt', '-ascii');
% line5 = load('Line5.txt', '-ascii');
% line6 = load('Line6.txt', '-ascii');
% line7 = load('Line7.txt', '-ascii');
% line8 = load('Line8.txt', '-ascii');
% line9 = load('Line9.txt', '-ascii');
% line10 = load('Line10.txt', '-ascii');
% line11 = load('Line11.txt', '-ascii');


inpath = './snapshots/';


tdelay = 1/24; % 24 frames per second

% parameters for plots
width = 0.85;
height = 0.5;
ygap = 0.003;
left_cor_x = 0.1;
left_cor_y = 0.02;

filename = 'displacements.gif';
figure(1);

perc = 1.0;
zclipx = 3;
zclipz = 3;
for isnap = 1:1:nsnap
    nstep = isnap*10;
    %fprintf('%5.5d \n',k);
    %======================================================================
    subplot(2,1,2);
    snapname = strcat(inpath,sprintf('%5.5d',nstep),'uz.su');
    snapz = readsu(snapname);
    imagesc(x, z, -snapz, perc*[-zclipz, zclipz]);
    hold on;
    % plot(line1(:,1), line1(:,2), 'k', 'linewidth', 0.8);
    % plot(line2(:,1), line2(:,2), 'k', 'linewidth', 0.8);
    % plot(line3(:,1), line3(:,2), 'k', 'linewidth', 0.8);
    % plot(line4(:,1), line4(:,2), 'k', 'linewidth', 0.8);
    % plot(line5(:,1), line5(:,2), 'k', 'linewidth', 0.8);
    % plot(line6(:,1), line6(:,2), 'k', 'linewidth', 0.8);
    % plot(line7(:,1), line7(:,2), 'k', 'linewidth', 0.8);
    % plot(line8(:,1), line8(:,2), 'k', 'linewidth', 0.8);
    % plot(line9(:,1), line9(:,2), 'k', 'linewidth', 0.8);
    % plot(line10(:,1), line10(:,2), 'k', 'linewidth', 0.8);
    % plot(line11(:,1), line11(:,2), 'k', 'linewidth', 0.8);
    %plot box
    %plot([x0 x1 x1 x0 x0], [z0 z0 z1 z1 z0], 'color', 'k', 'linewidth', 0.8);
    %%plot Moho
    %plot([x0 x1], [zMoho, zMoho], 'color', 'k', 'linewidth', 0.8);
    %%plot abnormal
    %plot([xb0 xb1 xb1 xb0 xb0], [zb0 zb0 zb1 zb1 zb0], 'color', 'k', 'linewidth', 0.8);
    %%plot abnormal
    %plot([xb3 xb4 xb4 xb3 xb3], [zb3 zb3 zb4 zb4 zb3], 'color', 'k', 'linewidth', 0.8);
    axis equal;
    axis([x0 x1 z0 z1]);
    set(gca, 'fontsize', 18, 'Fontname', 'Times New Roman');
    xlabel('Distance (km)', 'fontsize', 24);
    ylabel('Depth (km)', 'fontsize', 24);
    hold off;
    %str = sprintf('uz snapshot [%5.2f s]', nstep*dt);
    str = sprintf('uz [%5.2f s]', nstep*dt);
    title(str, 'fontsize', 32, 'Fontname', 'Times New Roman');
    set(gca, 'position', [left_cor_x left_cor_y width height]);
    %set(gcf,'position',[100 45 1200 850]);
    % colormap(gca, 'gray');
    % colormap(gca, 'hsv');
    % colormap(gca, 'jet');
    colormap(gca, cmap);
    %======================================================================
    %subplot(2,1,1, 'position', [left_cor_x left_cor_y+height width height]);
    subplot(2,1,1);
    snapname = strcat(inpath,sprintf('%5.5d', nstep),'ux.su');
    snapx = readsu(snapname);
    imagesc(x, z, snapx, perc*[-zclipx, zclipx]);
    hold on;
    % plot(line1(:,1), line1(:,2), 'k', 'linewidth', 0.8);
    % plot(line2(:,1), line2(:,2), 'k', 'linewidth', 0.8);
    % plot(line3(:,1), line3(:,2), 'k', 'linewidth', 0.8);
    % plot(line4(:,1), line4(:,2), 'k', 'linewidth', 0.8);
    % plot(line5(:,1), line5(:,2), 'k', 'linewidth', 0.8);
    % plot(line6(:,1), line6(:,2), 'k', 'linewidth', 0.8);
    % plot(line7(:,1), line7(:,2), 'k', 'linewidth', 0.8);
    % plot(line8(:,1), line8(:,2), 'k', 'linewidth', 0.8);
    % plot(line9(:,1), line9(:,2), 'k', 'linewidth', 0.8);
    % plot(line10(:,1), line10(:,2), 'k', 'linewidth', 0.8);
    % plot(line11(:,1), line11(:,2), 'k', 'linewidth', 0.8);
    %plot box
    %plot([x0 x1 x1 x0 x0], [z0 z0 z1 z1 z0], 'color', 'k', 'linewidth', 0.8);
    %%plot Moho
    %plot([x0 x1], [zMoho, zMoho], 'color', 'k', 'linewidth', 0.8);
    %%plot abnormal
    %plot([xb0 xb1 xb1 xb0 xb0], [zb0 zb0 zb1 zb1 zb0], 'color', 'k', 'linewidth', 0.8);
    %%plot abnormal
    %plot([xb3 xb4 xb4 xb3 xb3], [zb3 zb3 zb4 zb4 zb3], 'color', 'k', 'linewidth', 0.8);
    axis equal;
    axis([x0 x1 z0 z1]);
    set(gca, 'fontsize', 18, 'Fontname', 'Times New Roman');
    %set(gcf,'position',[400 50 1200 900]);
    xlabel('Distance (km) ', 'fontsize', 24);
    ylabel('Depth (km)', 'fontsize', 24);
    hold off;
    %str = sprintf('ux snapshot [%5.2f s]', nstep*dt);
    str = sprintf('ux [%5.2f s]', nstep*dt);
    title(str, 'fontsize', 32, 'Fontname', 'Times New Roman');
    set(gca, 'position', [left_cor_x left_cor_y+height+ygap width height]);
    % colormap(gca, 'gray');
    % colormap(gca, 'hsv');
    % colormap(gca, 'jet');
    colormap(gca, cmap);
    %======================================================================
    %set(gcf,'position',[100 45 1200 850]);
    set(gcf,'position',[100 45 1200 800]);
    %set(gcf,'position',[100 45 1200 1000]);
    set(gcf,'color','w')
    %======================================================================
    %pause(1/24); % 24 frames per second
    %pause(1/48); % 48 frames per second
    f = getframe(gcf);  
    imind = frame2im(f);
    [imind,cm] = rgb2ind(imind,256);
    if (isnap == 1)
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', tdelay);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', tdelay);
    end
end
