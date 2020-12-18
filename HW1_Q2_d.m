%Q2(d)

spectral_data=readmatrix('Problem 1.1 (Spectral data from the sun) - Sheet1.csv');
wavelength=spectral_data(:,1);
emission=spectral_data(:,2);
I=wavelength.*emission;
I_net=sum(I);
P=I_net/(3*10^8)