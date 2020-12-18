%EMA601 PS3
%MATLAB Code for numerical analysis of schrodinger equation

clear all

l=2; %angular wave number
n=3; %target radial wave number

n_range=linspace(n-0.02,n+0.02,100); %effective radial wave number

Et_range= -1./(n_range.^2); %E tilde range, corresponding to each radial wave number

Energy_range = 13.6.*Et_range; %actual energy range

dx=0.00001;
x_range=dx:dx:20; %x starts from dx and goes up to 100 and x=r/a_o

R_range = zeros(length(Et_range), length(x_range));
R_range(:,1)=1;
R_range(:,2)=1;

for i=1:1:length(Et_range)
    for j=1:1:(length(x_range)-2)
        R_range(i,j+2) = (2*R_range(i,j+1) + (l*(l+1)/x_range(j+1)^2-2/x_range(j+1)-Et_range(i))*R_range(i,j+1)*dx^2 - (1-dx/x_range(j+1))*R_range(i,j)) / (1+dx/x_range(j+1));
    end
end

%function to plot all R(x) values over the n* range
%clf(figure(1),'reset')
%for k=1:1:length(Et_range)
%    figure(1)
%    hold on
%    plot(x_range,R_range(k,:),'DisplayName',append('n^*=',num2str(n_range(k))),'Linewidth',1)
%end
%legend('Location', 'southeastoutside')
%title('R(x) vs x for various values of E')
%xlabel('x')
%ylabel('R(x)')
%ax = gca;
%ax.FontSize = 25;

diff=1e10;
E_choice=0;

for k=1:1:length(Et_range)
    if abs(R_range(k,end)) <= diff
        diff = abs(R_range(k,end));
        E_choice=k;
    end
end

n_range(E_choice) %effective n for the chosen solution
Energy_range(E_choice) %actual eigenenergy for the chosen solution

%function to plot the accepted R(x) solution
clf(figure(2),'reset')
figure(2)
plot(x_range,R_range(E_choice,:),'DisplayName',append('n^*=',num2str(n_range(E_choice))),'Linewidth',2)
legend('Location', 'southeastoutside')
title('R(x) vs x for accepted solution')
xlabel('x')
ylabel('R(x)')
ax = gca;
ax.FontSize = 25;

%function to plot the accepted rR(x) solution
r_range = 5.29177e-11.*x_range;
rR_range = r_range.*R_range(E_choice,:);
clf(figure(3),'reset')
figure(3)
plot(x_range,rR_range,'DisplayName',append('n^*=',num2str(n_range(E_choice))),'Linewidth',2)
legend('Location', 'southeastoutside')
title('r(x)R(x) vs x for accepted solution')
xlabel('x')
ylabel('r(x)R(x)')
ax = gca;
ax.FontSize = 25;

