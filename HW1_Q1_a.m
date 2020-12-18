%Script to find the wavelength at which intensity of blackbody radiation
%maximizes

%Defining parameters
c=2.99792*10^8; %Speed of light
h=6.626*10^(-34); %Planck's constant
kb=1.38*10^(-23); %Boltzmann constant

syms x 

sol_x = vpasolve((x/(1-exp(-x)))-5,x)

syms T %temperature
lambda=(h*c)/(sol_x*kb*T)

