load('Xm20to20um.mat')

%% useful constants
% permeability in H/um
mu_0 = 4 * pi * 1e-13;
 
% flux quanta in Wb 
Phi_0 = 2.068 * 1e-15;

xRange = x(end);
yRange = y(end);

dx = x(2) - x(1);
dy = y(2) - y(1);

%% zoom in to center part
dxp = .1;
dyp = .1;
% plot range
xP = 20;
yP = 20;
% plot range
xP = 20;
yP = 20;

% number of points to be plotted
Nx = round(xP/xRange * length(x)/2);
Ny = round(yP/yRange * length(x)/2);
% half of every vector
halfX = length(x)/2;
halfY = length(y)/2;
% range to interpolate to
xq = x(halfX - Nx):dxp:x(halfX + Nx + 2);
yq = y(halfX - Ny):dyp:y(halfY + Ny + 2);
[Xq,Yq] = meshgrid(xq,yq);
kx_useful = kx(halfX - Nx: halfX + Nx + 2);
ky_useful = ky(halfY - Ny: halfY + Ny + 2);

x_useful = x(halfX - Nx:halfX + Nx + 2);
y_useful = y(halfY - Ny:halfY + Ny + 2);

[X_useful, Y_useful] = meshgrid(x_useful, y_useful);
xq = xq - x(halfX + 1) * ones(1, length(xq));
yq = yq - y(halfY + 1) * ones(1, length(yq));

%% obtain expression in real space
for i = 1:length(x0)

hzk0 = hzk(:,:,i);

hz0 = ifft2(ifftshift(hzk0));
hz0 = fftshift((hz0))/(dx*dy);

hzk_useful = real(hzk(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY + Ny +2));

real_hz = real(hz0);
real_hz_useful = real_hz(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY+ Ny + 2);

real_hz_useful_interp = interp2(X_useful,Y_useful,real_hz_useful,Xq,Yq);
% convert to force
I = 3e-3; % A
XForce = 2 * Phi_0 * real_hz_useful_interp * 1e6 * I * 1e15; %fN

figure(i)
imagesc(xq, yq, XForce)
colormap jet
%colorbar 
c = colorbar;
c.Label.String = 'fN';
set(gca, 'FontSize',20)
%xlabel('x (\mum)')
%ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$\delta F_{x}$','Interpreter','latex')
hold on
%caxis([-8e-3 8e-3])
caxis([-1.2e-3 1.2e-3])
circle(x0(i), 0, 7.13)
circle(x0(i), 0, 3.57)
set(gca,'xtick',[])
set(gca,'ytick',[])
hold off
end

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,'k-','LineWidth', 3);
hold off
end