%% define the k-space
% real space range in um
dx = 1;
dy = 1;
xRange = 1500 + dx/2;
yRange = 1500 + dy/2;
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
[Kx, Ky] = meshgrid(kx,ky);
K = sqrt(Kx.^2 + Ky.^2);

%% define SQUID height and geometry
a = 7.13;
b = 3.57;
z_0 = 4.2;
alpha = 1;

%% center of SQUID
x0 = 5;
y0 = 2;
%% Pearl lenth and flux quanta
lambda_0 = 420;
Phi_0 = 2.068e-15; % N*m/A
%% get fields in k-space
% really hx(k)/I

%hxk
hxk = 1i * pi * a * Kx .* exp(-K * z_0 - 1i * Kx * x0 - 1i * Ky * y0) .* besselj(1, K * a) ./ (K .* (1 + lambda_0 * K));

%hyk
hyk = 1i * pi * a * Ky .* exp(-K * z_0- 1i * Kx * x0 - 1i * Ky * y0) .* besselj(1, K * a) ./ (K .* (1 + lambda_0 * K));

% only imag parts because real parts are zero
figure(288)
imagesc(kx,ky,imag(hxk))
colormap jet
c = colorbar;
c.Label.String = '\mum';
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_x^r(\vec{k})/I$','Interpreter','latex')

figure(287)
imagesc(kx,ky,imag(hyk))
colormap jet
c = colorbar;
c.Label.String = '\mum';
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gca, 'ydir', 'normal')
title('$h_y^r(\vec{k})/I$','Interpreter','latex')

%% Fourier transform for x and y 

hx0 = ifft2(ifftshift(hxk));
hx0 = fftshift((hx0))/(dx*dy);

hy0 = ifft2(ifftshift(hyk));
hy0 = fftshift((hy0))/(dx*dy);

figure(343)
imagesc(x,y,real(hx0))
c = colorbar;
c.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_x^r(\vec{r})/I$','Interpreter','latex')

figure(243)
imagesc(x,y,real(hy0))
colormap jet
c = colorbar;
c.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_y^r(\vec{r})/I$','Interpreter','latex')


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
hxk_useful = real(hxk(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY + Ny +2));
hyk_useful = real(hyk(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY + Ny +2));

figure(345)
imagesc(kx_useful, ky_useful, hxk_useful)
colormap jet
%colorbar
c = colorbar;
c.Label.String = '\mum';
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_x^r(\vec{k})/I$','Interpreter','latex')

figure(245)
imagesc(kx_useful, ky_useful, hyk_useful)
colormap jet
%colorbar
c = colorbar;
c.Label.String = '\mum';
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_y^r(\vec{k})/I$','Interpreter','latex')

x_useful = x(halfX - Nx:halfX + Nx + 2);
y_useful = y(halfY - Ny:halfY + Ny + 2);
[X_useful, Y_useful] = meshgrid(x_useful, y_useful);

real_hx = real(hx0);
real_hx_useful = real_hx(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY+ Ny + 2);

real_hy = real(hy0);
real_hy_useful = real_hy(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY+ Ny + 2);

figure(346)
imagesc(x_useful, y_useful, real_hx_useful)
colormap jet
%colorbar 
c = colorbar;
c.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_z^r(\vec{r})/I$','Interpreter','latex')
%title('real part, unit 1/um')

figure(446)
imagesc(x_useful, y_useful, real_hy_useful)
colormap jet
%colorbar 
c = colorbar;
c.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_y^r(\vec{r})/I$','Interpreter','latex')

%% interpolate center data
real_hx_useful_interp = interp2(X_useful,Y_useful,real_hx_useful,Xq,Yq);
real_hy_useful_interp = interp2(X_useful,Y_useful,real_hy_useful,Xq,Yq);

xq = xq - x(halfX + 1) * ones(1, length(xq));
yq = yq - y(halfY + 1) * ones(1, length(yq));

[Xq, Yq] = meshgrid(xq, yq);

figure(347)
imagesc(xq, yq, real_hx_useful_interp)
colormap jet
%colorbar 
c = colorbar;
c.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_x^r(\vec{r})/I$','Interpreter','latex')
hold on
circle(x0, y0, a);
circle(x0, y0, b);
hold off


figure(447)
imagesc(xq, yq, real_hy_useful_interp)
colormap jet
%colorbar 
c = colorbar;
c.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_y^r(\vec{r})/I$','Interpreter','latex')
hold on 
circle(x0, y0, a);
circle(x0, y0, b);
hold off


totalField = sqrt(real_hx_useful_interp.^2 + real_hy_useful_interp.^2);
figure(448)
imagesc(xq, yq, totalField)
colormap summer
c = colorbar;
c.Label.String = '1/\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'ydir', 'normal')
title('$h_{in}^r(\vec{r})/I$','Interpreter','latex')
hold on
quiver(Xq(5:25:end,5:25:end),Yq(5:25:end,5:25:end),real_hx_useful_interp(5:25:end,5:25:end),real_hy_useful_interp(5:25:end,5:25:end),.5,'r', 'LineWidth', 2)
circle(x0, y0, a);
circle(x0, y0, b);
hold off

I = 3e-3; % A
totalForce = 2 * Phi_0 * totalField * 1e6 * I * 1e15; %pF
figure(449)
imagesc(xq, yq, totalForce)
colormap summer
c = colorbar;
c.Label.String = 'pN';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
title('$F_{in}$','Interpreter','latex')
set(gca, 'ydir', 'normal')
hold on
quiver(Xq(1:25:end,1:25:end),Yq(1:25:end,1:25:end),real_hx_useful_interp(1:25:end,1:25:end),real_hy_useful_interp(1:25:end,1:25:end),.5,'r', 'LineWidth', 2)
circle(x0, y0, a);
circle(x0, y0, b);
hold off

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,'k-','LineWidth', 3);
hold off
end