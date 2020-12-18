%EMA601 PS3
%MATLAB Code for Q2d

clear all

dx=0.001;
x_range=dx:dx:1;

Z=11; %for sodium
r_c=1e-10; %r_c=0.1 nm for sodium


epo=8.845e-12; %vacuum permittivity
hbar=1.054e-34; %h/2pi
m=9.109e-31; %mass of electron
e=1.602e-19; %charge of electron
ao = (4*pi*epo*hbar^2)/(m*e^2);

n=4; %target radial wave number
l=1;

n_range=linspace(n-0.02,n+0.02,100); %effective radial wave number
Et_range= -1./(n_range.^2); %E tilde range, corresponding to each radial wave number
Energy_range = 13.6.*Et_range; %actual energy range

R_range = zeros(length(Et_range), length(x_range));
R_range(:,1)=1;
R_range(:,2)=1;

r_range=ao*x_range;

V_range = zeros(1,length(x_range));
Vr_range=zeros(1,length(r_range));

for i=1:1:(length(x_range))
    if r_range(i)<=r_c
        V_range(i)= -(2*Z)/(x_range(i)) + (2*ao*(Z-1))/r_c + (l*(l+1))/(x_range(i)^2);
        Vr_range(i) = -Z*e^2/(4*pi*epo*r_range(i)) + e^2*(Z-1)/(4*pi*epo*r_c);
    else
        V_range(i)= -2/x_range(i) + l*(l+1)/x_range(i)^2;
        Vr_range(i) = -e^2/(4*pi*epo*r_range(i));
    end
end

for i=1:1:length(Et_range)
    for j=1:1:(length(x_range)-2)
        R_range(i,j+2) = (2*R_range(i,j+1) + (V_range(j+1)-Et_range(i))*R_range(i,j+1)*dx^2 - (1-dx/x_range(j+1))*R_range(i,j)) / (1+dx/x_range(j+1));
    end
end

diff=1e10;
E_choice=0;

for k=1:1:length(Et_range)
    if abs(R_range(k,end)) <= diff
        diff = abs(R_range(k,end));
        E_choice=k;
    end
end

n_eff= n_range(E_choice) %effective n for the chosen solution
eigenenergy = Energy_range(E_choice) %actual eigenenergy for the chosen solution
defect = n_eff - n
