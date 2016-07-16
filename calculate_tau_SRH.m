%% Process sample lifetime data
clear all; close all; 
start = 'C:\Users\Mallory\Documents\Non-contact crucible\Australia experiments\June 22 2016\20-3\';

load([start 'averageTau.mat']); 
% data = dataSave{i}; 
% deltan = data(:,1); 
% tau = data(:,2); 
% 
% %Interpolate the lifetime data to use with the SRV later
% deltanq = logspace(13,17,500); 
% tauq = interp1(deltan,tau,deltanq); 

deltanq = deltanq;
tauq = tau_mean;

%% Process with surface and Richter

load('C:\Users\Mallory\Documents\Non-contact crucible\Australia experiments\June 22 2016\16-4-28-N\SRV_deltan.mat');

%Interpolate the SRV so that it matches the measured lifetime
SRVq = interp1(deltan,SRV,deltanq); 

D=11.84;
T = 300; 
N_dop = 2.8e15; %cm-3
type = 'n'; 
W = 0.0167; %cm

tau_surf =(W./(2.*SRVq))+((1/D).*((W/pi)^2)); %cm/s
tau_surf = tau_surf';

for i = 1:length(deltanq)
    tau_intr(i,1) = Richter(T,deltanq(i),N_dop,type);
end

figure;
loglog(deltanq,tauq.*1e6);
hold all; 
loglog(deltanq,tau_surf.*1e6); 
hold all; 
loglog(deltanq,tau_intr.*1e6); 

xlabel('Excess carrier density (cm^-^3)','FontSize',30);
ylabel('Lifetime (\mus)','FontSize',30);
axis([5e13 1e17 0 25000]);

tau_SRH = ((1./tau_mean)-(1./tau_intr)-(1./tau_surf)).^(-1);
hold all;
loglog(deltanq,tau_SRH.*1e6);
legend('Measured','Surface','Intrinsic','SRH');

save([start 'lifetime_breakdown.mat'],'tau_SRH','tau_surf','tauq','deltanq','tau_intr','N_dop','type','W');

%% Linearize the SRH term

%Plot the Murphy linearization for n-type
%Get sample parameters at specified temperature
[Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,N_dop,type); 
%Normalized carrier density
if type == 'p'
    X = (n0+deltanq)./(p0+deltanq);
elseif type == 'n'
    X = (p0+deltanq)./(n0+deltanq);
end
[m,n] = size(tau_SRH);
for i = 1:m
    figure;
    plot(X,tau_SRH(i,:).*1e6,'LineWidth',4);
    xlabel('X = p/n','FontSize',30);
    ylabel('SRH \tau (\mus)','FontSize',30);
    axis([0 1 0 2000]);
end
save('C:\Users\Mallory\Documents\Australia\Quasi mono NTNU\QSSPC\A44\linSRH.mat','tau_SRH','X','tau_surf','tauq','deltanq','tau_intr','N_dop','type','W');
% save('linSRH_SRVzero.mat','tau_SRH','X','tauq','deltanq','tau_intr','N_dop','type','W');
xlswrite('C:\Users\Mallory\Documents\Australia\Quasi mono NTNU\QSSPC\A44\linSRH.xlsx',[X,tau_SRH']);

%% Try fitting the linear lines 

clear all; close all;
load('FC94-23-4_summary.mat');
threshold = 1e-20; 

%Get rid of any NaN's because this will mess up the fitting
ix = find(isnan(tau_SRH)==1);
X(ix)=[];
tau_SRH(ix)=[];

[p,order,MSE,MSE_store]=fit_murphy(X,tau_SRH',threshold);