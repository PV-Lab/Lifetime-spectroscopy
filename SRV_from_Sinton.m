%%Load the data
clear all; close all; 
dirname = 'C:\Users\Mallory\Documents\Lifetime spectroscopy\Experiments\NOC-Si\July 17 2015\Sinton\15-7-14-N\Excel files'
process_xls_data(dirname);

%% Given the measured lifetime, determine the SRV
load('all_XLS_data.mat'); 

N_dop = 1.5e15; %cm-3
W = 0.0280; %cm
T = 300; %K
type = 'n';

%Get the intrinsic lifetime
tau_intr = zeros(length(deltanq),1);
%Calculate Richter model for the different doping levels
for i = 1:length(deltanq)
    tau_intr(i,1) = Richter(T,deltanq(i),N_dop,type);
end

%Get the surface lifetime
tau_surf = ((1./tauq_revised)-(1./tau_intr)).^(-1);
figure;
loglog(deltanq,tau_surf);
hold all;
loglog(deltanq,tau_intr);
hold all;
loglog(deltanq,tauq_revised);
legend('Surface','Intrinsic','Measured');

D = 11.97; 

SRV = W./((tau_surf-((1/D)*((W/pi)^2))).*2);

figure;
semilogx(deltanq,SRV); 
xlabel('Excess carrier density [cm^{-3}]');
ylabel('SRV [cm/s]');

save('SRV_deltan.mat','deltanq','SRV');

