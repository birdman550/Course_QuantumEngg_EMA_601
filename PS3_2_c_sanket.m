%EMA601 PS3
%MATLAB Code for Q2c

clear all

dx=0.001;
x_range=dx:dx:0.2;

Z=11; %for sodium
r_c=1e-10; %r_c=0.1 nm for sodium

l=2;
epo=8.845e-12; %vacuum permittivity
hbar=1.054e-34; %h/2pi
m=9.109e-31; %mass of electron
e=1.602e-19; %charge of electron
ao = (4*pi*epo*hbar^2)/(m*e^2);

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

%function V(x)
clf(figure(8),'reset')
figure(8)
hold on
plot(x_range,V_range,'DisplayName','V(x)','Linewidth',2)
plot(x_range,10^18*Vr_range,'DisplayName','V(r)*10^{18}','Linewidth',2)
legend('Location', 'southeastoutside')
title(append('V(x), V(r) vs x for Na at l=',num2str(l)))
xlabel('x')
ylabel('V(x), V(r)')
ax = gca;
ax.FontSize = 25;