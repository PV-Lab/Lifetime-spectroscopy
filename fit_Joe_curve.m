function [Joe,Joe_select] = fit_Joe_curve(deltan,tau,doping,T,type,fit_range,thickness)
%This function takes as input a lifetime and excess carrier density curve
%and does the fitting for the Joe value. 
%deltan = vector of injection levels [cm^-3]
%tau = vector of lifetimes [s], same size as deltan
%doping = sample doping level, measured at room temperature [cm^-3]
%T = sample temperature [K]
%type = 'p' or 'n'
%fit_range = range of injection levels to perform linear fit, expressed as
%a decimal percent (Sinton uses +/-30%, so input 0.30)
%Thickness = sample thickness, in cm

%First try the bulk lifetime with Auger
figure; 
plot(deltan,1./tau,'o'); 
%Then subtract Auger and plot it on the same plot
tau_intr = zeros(size(deltan));
for i = 1:length(deltan)
    tau_intr(i) = Richter(T,deltan(i),doping,type);
end
tau_bulk = 1./((1./tau)-(1./tau_intr)); 
hold all;
plot(deltan,1./tau_bulk,'o'); 

%Perform the fitting. 
minMCD = min(deltan); 
maxMCD = max(deltan); 
for i = 1:length(deltan); 
    minMCD_fit = deltan(i)*(1-fit_range);
    maxMCD_fit = deltan(i)*(1+fit_range); 
    indices = find(deltan>minMCD_fit & deltan<maxMCD_fit);
    deltan_fit = deltan(indices);
    tau_fit = 1./tau_bulk(indices); 
    p = polyfit(deltan_fit,tau_fit,1); 
    slope(i) = p(1); 
    intercept(i) = p(2); 
end

%Boltzmann constant
k_B = 8.61733238e-5; %eV/K  
%Bandgap in silicon
if T <= 303 & T >= 297
    Eg = 1.1242; %eV
    [NC,NV] = DOS_std(T); %cm^-3  
else
    [Eg] = Sze(T); %eV
    %Density of states effective mass model for consistency
    [NC,NV] = DOS_em(T); %cm^-3  
end
ni2 = NC*NV*exp(-Eg/(k_B*T)); %cm^-6
q = 1.602e-19; %C
Joe = slope.*q.*ni2.*thickness./2;
%We want to be one order of magnitude above the doping concentration
if doping*10>max(deltan)
    %Then we choose the fit at highest injection
    for_average = Joe(1:5); 
    Joe_select = mean(for_average); 
    slope_select = mean(slope(1:5));
    intercept_select = mean(intercept(1:5));
%     Joe_select= Joe(find(deltan==max(deltan)));
%     slope_select = slope(find(deltan==max(deltan)));
%     intercept_select = intercept(find(deltan==max(deltan)));
else
    %We choose the value close to one order of magnitude higher than
    %doping
%     diff = abs(doping*10 - doping); 
    for_average = Joe(find(deltan>doping*10));
    Joe_select = mean(for_average); 
    slope_select = mean(slope(find(deltan>doping*10)));
    intercept_select = mean(intercept(find(deltan>doping*10)));
%     Joe_select = Joe(find(diff==min(diff))); 
%     slope_select = slope(find(diff==min(diff)));
%     intercept_select = intercept(find(diff==min(diff)));
end 
%Plot the chosen fit
hold all; 
plot(deltan,slope_select.*deltan+intercept_select,'k-','LineWidth',3);
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Inverse lifetime [s^-^1]','FontSize',20);
legend('Auger not removed','Auger removed','Chosen fit'); 
disp('Fitting will be performed on data with Auger removed.');


figure;
loglog(deltan,Joe,'o'); 
hold all; 
loglog([doping doping],[-1 1],'k-');
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('J_o_e [A/cm^2]','FontSize',20);
title('Possible J_o_e values','FontSize',30); 


