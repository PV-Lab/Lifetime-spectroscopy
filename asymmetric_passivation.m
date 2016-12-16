%{
MIT License

Copyright (c) [2016] [Mallory Ann Jensen, jensenma@alum.mit.edu]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%}

%%Plot lifetime for different versions of asymmetrically passivated wafers
clear all; close all; 

%%First for a sample with symmetrically passivated emitters
Joe = 5e-14; %A/cm2
type = 'p'; 
doping = 1.5e16; %/cm3
deltan = logspace(13,18,500); %/cm3
q = 1.602e-19; %C
W = .0180; %cm, wafer thickness
T = [300 350 400 450]; %K
%Get the intrinsic carrier concentration for each temperature
%Call standard density of states model
%Boltzmann constant
k_B = 8.61733238e-5; %eV/K
ni2 = zeros(size(T)); 
for tempi = 1:length(T)
    [NC,NV] = DOS_std(T(tempi)); %cm^-3  
    %Bandgap in silicon
    [Eg] = Sze(T(tempi)); %eV
    ni2(tempi) = NC*NV*exp(-Eg/(k_B*T(tempi))); %cm^-6
end
%Assume that Joe remains constant with temperature and calculate surface
%lifetime
tau_surf = zeros(length(deltan),length(T)); 
figure;
for tempi = 1:length(T)
    tau_surf_hold = (2.*Joe.*(doping+deltan)./(q*ni2(tempi)*W)).^(-1); 
    tau_surf(:,tempi) = tau_surf_hold'; 
    h(tempi) = loglog(deltan,tau_surf_hold,'LineWidth',3);
    hold all; 
    %Also get Richter lifetime
    for inji = 1:length(deltan)
        tau_intr(inji) = Richter(T(tempi),deltan(inji),doping,type);
    end
    hold all; 
    loglog(deltan,tau_intr,'k--');
end

xlabel('Excess carrier density [cm^-^3]','FontSize',30); 
ylabel('Lifetime [s]','FontSize',30);
legend(h,num2str(T'));
title('Symmetric passivated emitter','FontSize',30);

%%Now a sample with asymmetric passivation
Joe = 5e-14; %A/cm2
type = 'p'; 
doping = 1.5e16; %/cm3
deltan = logspace(13,18,500); %/cm3
q = 1.602e-19; %C
W = .0180; %cm, wafer thickness
T = [300 350 400 450]; %K
S = 20; %cm/s
%Get the intrinsic carrier concentration for each temperature
%Call standard density of states model
%Boltzmann constant
k_B = 8.61733238e-5; %eV/K
ni2 = zeros(size(T)); 
for tempi = 1:length(T)
    [NC,NV] = DOS_std(T(tempi)); %cm^-3  
    %Bandgap in silicon
    [Eg] = Sze(T(tempi)); %eV
    ni2(tempi) = NC*NV*exp(-Eg/(k_B*T(tempi))); %cm^-6
end
%Assume that Joe remains constant with temperature and calculate surface
%lifetime
tau_surf = zeros(length(deltan),length(T)); 
figure;
for tempi = 1:length(T)
    tau_surf_hold = ((Joe.*(doping+deltan)./(q*ni2(tempi)*W))+(S/W)).^(-1); 
    tau_surf(:,tempi) = tau_surf_hold'; 
    h(tempi) = loglog(deltan,tau_surf_hold,'LineWidth',3);
    hold all; 
    %Also get Richter lifetime
    for inji = 1:length(deltan)
        tau_intr(inji) = Richter(T(tempi),deltan(inji),doping,type);
    end
    hold all; 
    loglog(deltan,tau_intr,'k--');
end

xlabel('Excess carrier density [cm^-^3]','FontSize',30); 
ylabel('Lifetime [s]','FontSize',30);
legend(h,num2str(T'));
title('Asymmetric passivated emitter','FontSize',30);

