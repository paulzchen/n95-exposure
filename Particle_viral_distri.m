%This script calculates the expected percentage of airborne particles 
%containing SARS-CoV-2 or A(H1N1)pdm09 virions for the relevant size range
%(0.1 - 5 um) according to the maximum probability case, as described in 
%the manuscript. Also, we apply a conservative estimate and consider each
%copy in the clinical viral loads to represent a viable virion. 

%This script also calculates virion partioning upon airborne particle 
%generation (based on the Poisson distribution) for SARS-CoV-2 and 
%A(H1N1)pdm09 based on the mean viral loads.

clear all
close all
clc

%Define constants
p = 997e3; %Density of water at room temperature (g/m3)
Vml_Vg = 1; %Conversion factor (volume) from g to ml for water at room
            %temperature
C_SARS2 = 9.25e5; %Viral load (average) in respiratory specimens for 
               %SARS-CoV-2 (copies per ml)
C_H1N1 = 2.63e6; %Viral load (average) in respiratory specimens for 
                 %A(H1N1)pdm09 (copies per ml)
d = linspace(0.1,5,201)./0.44.* 1e-6; %Size range of hydrated airborne
                                      % particles during generation (m)
d_deh = d.*0.44*1e6; %Size range of dehydrated airborne particles (um)

V = p* Vml_Vg* 4/3* pi* (d./2).^3;
pSARS2 = C_SARS2.* V;
pH1N1 = C_H1N1.* V;


%Approximate the probability distribution of virion partioning
lamda = [pSARS2(1), pSARS2(49), pSARS2(end), pH1N1(1), pH1N1(end)];
v = [0 1 2 3 4];
Pv_SARS2_0_1um = zeros(length(v),1);
Pv_SARS2_5um = zeros(length(v),1);
Pv_SARS2_11_3um = zeros(length(v),1);
Pv_H1N1_0_1um = zeros(length(v),1);
Pv_H1N1_11_3um = zeros(length(v),1);

for k = 1 : length(v)
    Pv_SARS2_0_1um(k) = (lamda(1)^v(k))* exp(-lamda(1)) / ...
        factorial(v(k));
    Pv_SARS2_5um(k) = (lamda(2)^v(k))* exp(-lamda(2)) / ...
        factorial(v(k));
    Pv_SARS2_11_3um(k) = (lamda(3)^v(k))* exp(-lamda(3)) / ...
        factorial(v(k));
    Pv_H1N1_0_1um(k) = (lamda(4)^v(k))* exp(-lamda(4)) / ...
        factorial(v(k));
    Pv_H1N1_11_3um(k) = (lamda(5)^v(k))* exp(-lamda(5)) / ...
        factorial(v(k));
end


subplot(2,2,1)
plot(d_deh,V)
xlabel('Airborne particle size (um)')
ylabel('Volume (ml)')

subplot(2,2,2)
plot(d_deh,pSARS2.*100)
title('SARS-CoV-2')
xlabel('Airborne particle size (um)')
ylabel('Distribution (%)')

subplot(2,2,4)
plot(d_deh,pH1N1.*100)
title('A(H1N1)pdm09')
xlabel('Airborne particle size (um)')
ylabel('Distribution (%)')

figure;
subplot(2,3,1)
plot(v,Pv_SARS2_0_1um*100)
title('Probability that 0.1-um particle contains SARS-CoV-2')
xlabel('Number of virions')
ylabel('Distribution (%)')

subplot(2,3,2)
plot(v,Pv_SARS2_5um*100)
title('Probability that 5-um particle contains SARS-CoV-2')
xlabel('Number of virions')
ylabel('Distribution (%)')

subplot(2,3,3)
plot(v,Pv_SARS2_11_3um*100)
title('Probability that 11.3-um particle contains SARS-CoV-2')
xlabel('Number of virions')
ylabel('Distribution (%)')

subplot(2,3,5)
plot(v,Pv_H1N1_0_1um*100)
title('Probability that 0.1-um particle contains A(H1N1)pdm09')
xlabel('Number of virions')
ylabel('Distribution (%)')

subplot(2,3,6)
plot(v,Pv_H1N1_11_3um*100)
title('Probability that 11.3-um particle contains A(H1N1)pdm09')
xlabel('Number of virions')
ylabel('Distribution (%)')