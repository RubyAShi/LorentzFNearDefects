% define real space range in um
dx = .1;
dy = .1;
xRange =600 + dx/2;
yRange =600 + dy/2;
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

a = 7.13;
b = 3.57;
z_0 = 4.4;
Lambda_0 = 357;

[Kx, Ky] = meshgrid(kx,ky);
k = sqrt(Kx.^2 + Ky.^2);
% plot hkz/I
hzk = -pi * a * exp(-2 * k * z_0) .* besselj(1, k * a) ./ (1 + Lambda_0 * k);

plot_range = 10;
% number of points to be plotted
N = round(plot_range/xRange * length(x)/2);
% half of every vector
half = length(x)/2;

figure(288)
imagesc(kx(half - N: half + N),ky(half - N: half + N),hzk(half - N: half + N, (half - N: half + N)))
colormap jet
%colorbar
a = colorbar;
a.Label.String = '\mum';
set(gca, 'FontSize',20)
xlabel('kx (1/\mum)')
ylabel('ky (1/\mum)')
set(gcf,'Position',[100 100 500 427])
title('$h_z^r(\vec{k})/I$','Interpreter','latex')

hz = ifft2(ifftshift(hzk));
hz = fftshift((hz))/(dx*dy);
real_hz = real(hz);
figure(343)
imagesc(x(half - N: half + N),y(half - N: half + N),real_hz(half - N: half + N, half - N: half + N))
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

imag_hz = imag(hz);
figure(344)
imagesc(x(half - N: half + N),y(half - N: half + N),imag_hz(half - N: half + N, half - N: half + N))
colormap jet
%colorbar 
a = colorbar;
a.Label.String = '\mum';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
title('imag part')

PickupLoop = 0*X;
R = sqrt(X.^2 + Y.^2);
PickupLoop(R - b <.001)=1;
figure(356)
imagesc(x(half - N: half + N),y(half - N: half + N),PickupLoop(half - N: half + N, half - N: half + N))
shading flat;
axis image;
colorbar;
colormap bone;
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'FontSize',20)

% try a faster method
xPU = -4:.1:4;
yPU = -4:.1:4;
[XPU, YPU] = meshgrid(xPU, yPU);
RPU = sqrt(XPU.^2 + YPU.^2);
PU = 0 * XPU;
PU(RPU - b <.001)=1;
figure(456)
imagesc(xPU,yPU,PU)
shading flat;
axis image;
colorbar;
colormap bone;
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'FontSize',20)

% find the -4 to 4 range
PU_range = 4;
% number of points to be plotted
NPU = round(PU_range/xRange * length(x)/2);

hz_PU = real_hz(half - NPU: half + NPU, half - NPU: half + NPU);
fluxPU = hz_PU .*PU * dx * dy;

fluxIMG = real(hz).* PickupLoop * dx * dy;
figure(357)
imagesc(x(half - N: half + N),y(half - N: half + N),fluxIMG(half - N: half + N, half - N: half + N))
shading flat;
axis image;
%colorbar;
a = colorbar;
a.Label.String = '\mum';
colormap jet;
%title('unit um')
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'FontSize',20)
title('$\Phi/I$','Interpreter','latex')

% unit um
%flux_um = sum(sum(fluxIMG));
%save('trial.mat','hzk', 'kx','ky', 'x', 'y', 'hz');

flux_um = sum(sum(fluxPU));
% permeability in H/um
mu_0 = 4 * pi * 1e-13;

% flux quanta in Wb 
Phi_0 = 2.068 * 1e-15;

% flux in Phi_0/A
flux = flux_um * mu_0/Phi_0;
