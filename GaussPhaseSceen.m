clc
close all
clear all
%Система координат
X=[-200:1:200]; %метры 1e-3
Y=[-200:1:200]; %метры 1e-3
[X Y]=meshgrid(X,Y);

[phi,ro] = cart2pol(X,Y);


w0=0.001; % beam spot m
lyam=532e-9; % длинна волны m; 
k=(2*pi)/lyam; %wavenumber; 


%ПЕРЕМЕННЫЕ НАКАЧКА, РАДИАЛЬНЫЕ И АЗИМУТАЛЬНЫЕ ИНДЕКСЫ

L1=2;
p1=0;
% Lp=0;
%%%%%%%%   по формуле 3 , при z=0
syms m
PolinomL=symsum((-1)^m.*((factorial(L1+p1))./((factorial(p1-m))*(factorial(L1+m))*(factorial(m))))*(2*ro.^2/w0^2).^m, m, 0, p1); %полином Лаггера формула (2)
P=(exp(-ro.^2/(w0^2))).*(exp(-1i*L1.*phi)); 
Lg=(sqrt(2*factorial(p1)/(pi*(factorial(L1+p1)))))*(1/w0)*((sqrt(2)*ro./w0).^L1).*P.*PolinomL; %

% s=class(Lg)
Lg=double(Lg); %Convert Symbolic Number to Double Precision
% C3=char(C)
% 
% mesh(X,Y,abs(Lg));
% set(gca,'view',[0 90]);grid on;

% % Параметры
I_0=1; %начальная интенсивноть 
g=100;

r0=0.0005; % радиус окружности в плоскости фокусировки
% sigma=0.005; %радиус окружности пучка

% Matlab code
sz=200; % size
% generate the power spectral density values
cx=(-sz:sz);
mx=(ones(2*sz+1,1)*cx).^2;
mr=sqrt(mx+transpose(mx));
psd=0.023*mr.^(-11/3);
psd(sz+1,sz+1)=0;
% generate the random numbers with Gaussian statistics
randomcoeffs=randn(2*sz+1)+1i*randn(2*sz+1);
% phase screen!
phasescreen=real(fft2(fftshift(sqrt(psd).*randomcoeffs)));
figure
mesh(phasescreen)


% I=I_0*exp((-(X.^2+Y.^2))/(2*(g^2)));; %функция гаусс пучка
% figure
% plot(I)


% I1=phasescreen.*Lg;
% figure
% mesh(I1)
% mesh(I)