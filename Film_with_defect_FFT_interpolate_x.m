%% define the k-space
% real space range in um
dx = 1;
dy = 1;
xRange = 400 + dx/2;
yRange = 400 + dy/2;
x = -xRange:dx:xRange;
y = -yRange:dy:yRange;
[X,Y] = meshgrid(x,y);
% defind the k-space
% dkx = 2pi/(range of x)
dkx = pi/xRange;
dky = pi/yRange;
% (range of kx) = 2pi/dx
kxRange = pi/dx;
kyRange = pi/dy;
kx = linspace(-kxRange, kxRange - dkx, length(x));
ky = linspace(-kyRange, kyRange - dky, length(y));
% SQUID location 
x0 = 3.5;
y0 = 0;
[Kx, Ky, X0] = meshgrid(kx,ky,x0);

%% define SQUID geometry 
a = 7.13;
b = 3.57;

%% define defect strength
alpha = 1;

%% other useful constants
% permeability in H/um
mu_0 = 4 * pi * 1e-13;

% flux quanta in Wb 
Phi_0 = 2.068 * 1e-15;

%% get hx
% really hx(k)/I
hzk = alpha^2 * a * arrayfun(@(kx, ky, x0) integral(@(qx) M(kx,ky,qx,x0),-inf, inf, 'RelTol',1e-5), Kx, Ky, X0);
%hzk = 2 * 1i * Kx .* pi * a * exp(-K * z_0) .* besselj(1, K * a) ./ (K .* (1 + lambda_0 * K));

hzk0 = hzk(:,:,1);
figure(288)
imagesc(kx,ky,real(hzk0))
colormap jet
colorbar
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gca, 'ydir', 'normal')
set(gcf,'Position',[100 100 500 427])
title('unit um')

hz0 = ifft2(ifftshift(hzk0));
hz0 = fftshift((hz0))/(dx*dy);
figure(343)
imagesc(x,y,real(hz0))
colormap jet
colorbar 
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gca, 'ydir', 'normal')
set(gcf,'Position',[100 100 500 427])
title('real part, unit 1/um')

figure(243)
imagesc(x,y,real(hz0))
colormap jet
colorbar 
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_z^r(\vec{r})/I$','Interpreter','latex')

% figure(344)
% imagesc(x,y,imag(hz0))
% colormap jet
% colorbar 
% set(gca, 'FontSize',20)
% xlabel('x (\mum)')
% ylabel('y (\mum)')
% set(gcf,'Position',[100 100 500 427])
% title('imag part')

%% zoom in to center part
dxp = .1;
dyp = .1;
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
hzk_useful = real(hzk(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY + Ny +2));

figure(345)
imagesc(kx_useful, ky_useful, hzk_useful)
colormap jet
%colorbar
c = colorbar;
c.Label.String = '\mum';
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gca, 'ydir', 'normal')
set(gcf,'Position',[100 100 500 427])
title('$h_x^r(\vec{k})/I$','Interpreter','latex')

x_useful = x(halfX - Nx:halfX + Nx + 2);
y_useful = y(halfY - Ny:halfY + Ny + 2);
[X_useful, Y_useful] = meshgrid(x_useful, y_useful);
real_hz = real(hz0);
real_hz_useful = real_hz(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY+ Ny + 2);
figure(346)
imagesc(x_useful, y_useful, real_hz_useful)
colormap jet
%colorbar 
c = colorbar;
c.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gca, 'ydir', 'normal')
set(gcf,'Position',[100 100 500 427])
title('$h_x^r(\vec{r})/I$','Interpreter','latex')
hold on
circle(x0, y0, a);
circle(x0, y0, b)
hold off
%title('real part, unit 1/um')

%% interpolate center data
real_hz_useful_interp = interp2(X_useful,Y_useful,real_hz_useful,Xq,Yq);

xq = xq - x(halfX + 1) * ones(1, length(xq));
yq = yq - y(halfY + 1) * ones(1, length(yq));

figure(347)
imagesc(xq, yq, real_hz_useful_interp)
colormap jet
%colorbar 
c = colorbar;
c.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
hold on
circle(x0, y0, a);
circle(x0, y0, b);
title('$h_x^r(\vec{r})/I$','Interpreter','latex')
hold off

%% convert to force
I = 3e-3; % A
XForce = 2 * Phi_0 * real_hz_useful_interp * 1e6 * I * 1e15; %pF

figure(447)
imagesc(xq, yq, XForce)
colormap jet
%colorbar 
c = colorbar;
c.Label.String = 'fN';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
hold on
circle(x0, y0, a);
circle(x0, y0, b);
title('$F_{x}$','Interpreter','latex')
hold off
%% function 
function y = M(kx, ky, qx, x0)
lambda_0 = 420;
z_0 = 4.2;
a = 7.13;
Q = sqrt(qx.^2 + ky.^2);
k = sqrt(kx.^2 + ky.^2);
xPrime = 7;
y0 = 0;
% plot hkz/I
%y = exp(-(Q + k)* z_0 + 1i * (qx-kx) * x0) .* besselj(1, Q * a) .* (qx .* kx + ky.^2) ./ (2 * (1 + lambda_0 * Q) .* (1 + lambda_0 * k).* Q);
y = 1i .* kx .* exp(-Q * z_0 - 1i * qx * x0- 1i * ky * y0) .* besselj(1, Q * a) .* (qx .* kx + ky.^2) .* cos((kx - qx)*xPrime) ./ (k .* (1 + lambda_0 * Q) .* (1 + lambda_0 * k).* Q);
end

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,'k-','LineWidth', 3);
hold off
end