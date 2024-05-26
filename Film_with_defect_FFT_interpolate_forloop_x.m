%% defind the k-space
% define real space range in um
dx = 1;
dy = 1;
xRange = 800 + dx/2;
yRange = 800 + dy/2;
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
% (x0,y0) here is the center of the SQUID
% due to symmetry, only x0 is changed
x0 = -20:2:20;
y0 = 0;
[Kx, Ky, X0] = meshgrid(kx,ky,x0);

%% SQUID geometry 
a = 7.13;
b = 3.57;

%% defect strength and location
alpha = 1;
xPrime = 7;
%% other useful constants
% permeability in H/um
mu_0 = 4 * pi * 1e-13;
 
% flux quanta in Wb 
Phi_0 = 2.068 * 1e-15;

%% get Hx
% really hx(k)/I
hzk = alpha^2 * a * arrayfun(@(kx, ky, x0) integral(@(qx) M(kx,ky,qx,x0),-inf, inf, 'RelTol',1e-5), Kx, Ky, X0);

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

x_useful = x(halfX - Nx:halfX + Nx + 2);
y_useful = y(halfY - Ny:halfY + Ny + 2);

[X_useful, Y_useful] = meshgrid(x_useful, y_useful);
xq = xq - x(halfX + 1) * ones(1, length(xq));
yq = yq - y(halfY + 1) * ones(1, length(yq));

%% save simulation in k-space
save('Xm20to20um.mat','hzk', 'kx','ky', 'x', 'y', 'a', 'b', 'alpha', 'x0', 'y0', 'xPrime');

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
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$F_{x}$','Interpreter','latex')
end

%% function 
function y = M(kx, ky, qx, x0)
lambda_0 = 400;
z_0 = 4;
a = 7.13;
Q = sqrt(qx.^2 + ky.^2);
k = sqrt(kx.^2 + ky.^2);
xPrime = 7;
y0 = 0;
% plot hkz/I
%y = exp(-(Q + k)* z_0 + 1i * (qx-kx) * x0) .* besselj(1, Q * a) .* (qx .* kx + ky.^2) ./ (2 * (1 + lambda_0 * Q) .* (1 + lambda_0 * k).* Q);
y = 1i .* kx .* exp(-Q * z_0 - 1i * qx * x0- 1i * ky * y0) .* besselj(1, Q * a) .* (qx .* kx + ky.^2) .* cos((kx - qx)*xPrime) ./ (k .* (1 + lambda_0 * Q) .* (1 + lambda_0 * k).* Q);
end


