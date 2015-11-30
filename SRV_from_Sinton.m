%%Load the data
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\2015 Oxygen-State Study\NOC samples\NOC 17 21 22 study\November 2015\One Day after Anneal\FZ';
process_xls_data(dirname);
%dataSave is a cell structure with probably one entry. The first column is
%injection level. The second column is lifetime. 

%% Given the measured lifetime, determine the SRV
% load('all_XLS_data.mat'); 
% data = dataSave{1}; %this assumes you're only looking at one sample
% deltan = data(:,1);
% tau = data(:,2);

load('average_data.mat');
tau = tau_mean; 
deltan = deltanq;

N_dop = 1.7e15; %cm-3
W = 0.0280; %cm
T = 300; %K
type = 'n';

%Get the intrinsic lifetime
tau_intr = zeros(length(deltan),1);
%Calculate Richter model for the different doping levels
for i = 1:length(deltan)
    tau_intr(i,1) = Richter(T,deltan(i),N_dop,type);
end

%Get the surface lifetime
tau_surf = ((1./tau)-(1./tau_intr)).^(-1);
figure;
loglog(deltan,tau_surf);
hold all;
loglog(deltan,tau_intr);
hold all;
loglog(deltan,tau);
legend('Surface','Intrinsic','Measured');

D = 11.95; 

SRV = W./((tau_surf-((1/D)*((W/pi)^2))).*2);

figure;
semilogx(deltan,SRV); 
xlabel('Excess carrier density [cm^{-3}]');
ylabel('SRV [cm/s]');

save('SRV_deltan.mat','deltan','SRV');

