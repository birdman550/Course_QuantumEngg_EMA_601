%EMA601 Mid-term project
%MATLAB scriptt to simulate sodium (Na-11) atom BECs
clear all

h=6.626e-34; %Planck's const in SI units
hb = h/(2*pi); %Reduced Planck's const in SI units
c=3e06; %speed of light in SI units
kb =  1.38e-23; %Boltzmann const in SI units

%defining initial parameters
a = 2.9e-9; %scattering length of sodium
wx = 2*pi * 250; %x-axis: angular oscillation freq
wy = wx; %y-axis: angular oscillation freq
wz = 2*pi * 16; %z-axis: angular oscillation freq
aho = 2e-6; %harmonic oscillator length
N = 1e06; %target number of atoms in the condensate


%determining secondary parameters
w = (wx*wy*wz)^(1/3) %average angular oscillation freq
M = hb/(w*aho^2); %mass of a sodium atom

%determining properties
Tc = hb*w*(N^(1/3))/kb %critical temperature for the trap
nc = 2.6/((h/sqrt(2*pi*M*kb*Tc))^3) %critical density of the trap; at critical temperature
mu = hb*w*0.5*(15*N*a/aho)^(2/5) %chemical potential of the atoms
Rx = sqrt(2*mu/(M*(wx^2))) %x-axis: radial extent of the condensate
Ry = sqrt(2*mu/(M*(wy^2))); %y-axis: radial extent of the condensate
Rz = sqrt(2*mu/(M*(wz^2))) %z-axis: axial extent of the condensate
g = 4*pi*(hb^2)*N*a/M; %coupling const for inter-atomic interaction

%[X,Y,Z] = meshgrid(-100:1:100,-100:1:100,-100:1:100); %X,Y,Z defined from -100 to 100 um, step size of 1 um
[X,Y] = meshgrid(-Rx:1e-7:Rx,-Ry:1e-7:Ry);
Z=0;
no = N*mu/g %number density at the center
n_range = no.*(1 - (X.^2)./(Rx^2) - (Y.^2)./(Ry^2) - (Z.^2)./(Rz^2));
n_range(n_range<0)=0;

clf(figure(1),'reset')
figure(1)
surf(X,Y,n_range, 'edgecolor', 'none')
title('number density in the x-y plane')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Number density (n) [m^{-3}]')
ax = gca;
ax.FontSize = 30;
%size(n_range);