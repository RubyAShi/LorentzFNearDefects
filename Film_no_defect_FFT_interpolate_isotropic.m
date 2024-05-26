%% define the k-space
% real space range in um
dx = 1;
dy = 1;
xRange = 1000 + dx/2;
yRange = 1000 + dy/2;
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
x0 = 0;
[Kx, Ky, X0] = meshgrid(kx,ky,x0);
K = sqrt(Kx.^2 + Ky.^2);

%% define SQUID height and geometry
a = 7.13;
b = 3.57;
z_0 = 4.2;
alpha = 1;

%% Pearl lenth
lambda_0 = 420;

%% get fields in k-space
% really hx(k)/I

%hxk
%hzk = 1i * pi * a * Kx .* exp(-K * z_0) .* besselj(1, K * a) ./ (K .* (1 + lambda_0 * K));

%hyk
hzk = 1i  * pi * a * Ky .* exp(-K * z_0) .* besselj(1, K * a) ./ (K .* (1 + lambda_0 * K));

%hzk at sample surface
%hzk =  -pi * a * exp(-K * z_0) .* besselj(1, K * a) ./ (1 + lambda_0 * K);

title_str = '$h_y^r(\vec{r})/I$';

hzk0 = hzk(:,:,1);
figure(288)
imagesc(kx,ky,real(hzk0))
colormap jet
colorbar
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gcf,'Position',[100 100 500 427])
title('unit um real part')

figure(287)
imagesc(kx,ky,imag(hzk0))
colormap jet
colorbar
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gcf,'Position',[100 100 500 427])
title('unit um imag part')

hz0 = ifft2(ifftshift(hzk0));
hz0 = fftshift((hz0))/(dx*dy);
figure(343)
imagesc(x,y,real(hz0))
colormap jet
colorbar 
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
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
title('$h_z^r(\vec{r})/I$','Interpreter','latex')

figure(344)
imagesc(x,y,imag(hz0))
colormap jet
colorbar 
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
title('imag part')

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
a = colorbar;
a.Label.String = '\mum';
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gcf,'Position',[100 100 500 427])
title('$h_z^r(\vec{k})/I$','Interpreter','latex')

x_useful = x(halfX - Nx:halfX + Nx + 2);
y_useful = y(halfY - Ny:halfY + Ny + 2);
[X_useful, Y_useful] = meshgrid(x_useful, y_useful);
real_hz = real(hz0);
real_hz_useful = real_hz(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY+ Ny + 2);
figure(346)
imagesc(x_useful, y_useful, real_hz_useful)
colormap jet
%colorbar 
a = colorbar;
a.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
title('$h_z^r(\vec{r})/I$','Interpreter','latex')
%title('real part, unit 1/um')

%% interpolate center data
real_hz_useful_interp = interp2(X_useful,Y_useful,real_hz_useful,Xq,Yq);

xq = xq - x(halfX + 1) * ones(1, length(xq));
yq = yq - y(halfY + 1) * ones(1, length(yq));

figure(347)
imagesc(xq, yq, real_hz_useful_interp)
colormap jet
%colorbar 
a = colorbar;
a.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
title(title_str,'Interpreter','latex')

