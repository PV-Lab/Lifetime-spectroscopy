%% Process raw lifetime data
clear all; close all; 
deltanq = logspace(13,17,500); 

load('all_XLS_data_broken.mat'); 

together = figure;

for i = 1:length(fileListShort)
    data = dataSave{i};
    deltan = data(:,1); 
    tau = data(:,2); 
    
    deltan = flipud(deltan);
    tau = flipud(tau); 
    [deltan,IX] = sort(deltan); 
    tau = tau(IX); 
    store_index = []; 
    count = 1; 
    for j = 1:length(deltan)
        index = find(deltan==deltan(j)); 
        [m,n] = size(index); 
        if m>1
            store_index(count:count+m-2) = index(2:end); 
            count = count+1; 
        end
    end
    deltan(store_index) = [];
    tau(store_index) = [];
    
    loglog(deltan,tau.*1e6,'.','MarkerSize',5);
    hold all; 
    
    tauq_broken(:,i) = interp1(deltan,tau,deltanq);
end

tau_mean_broken = nanmean(tauq_broken,2); 

loglog(deltanq,tau_mean_broken.*1e6,'k-','LineWidth',2); 
hold all; 

xlabel('Excess carrier density (cm^-^3)','FontSize',30);
ylabel('Lifetime (\mus)','FontSize',30);
axis([5e13 1e17 0 8000]);

load('all_XLS_data_NOTbroken.mat'); 

for i = 2
    data = dataSave{i};
    deltan = data(:,1);
    tau = data(:,2);
    
    deltan = flipud(deltan);
    tau = flipud(tau); 
    [deltan,IX] = sort(deltan); 
    tau = tau(IX); 
    store_index = []; 
    count = 1; 
    for j = 1:length(deltan)
        index = find(deltan==deltan(j)); 
        [m,n] = size(index); 
        if m>1
            store_index(count:count+m-2) = index(2:end); 
            count = count+1; 
        end
    end
    deltan(store_index) = [];
    tau(store_index) = [];
    
    loglog(deltan,tau.*1e6,'.','MarkerSize',5);
    hold all;
    
    tauq(:,i) = interp1(deltan,tau,deltanq);
end

% tau_mean = nanmean(tauq,2); 
tau_mean = tauq(:,i);

loglog(deltanq,tau_mean.*1e6,'k-','LineWidth',2); 

%Find the multiplier between the two data sets
index = find(deltanq>=1e15);
tau_mean_fit = tau_mean(index);
tau_mean_broken_fit = tau_mean_broken(index);

%Get rid of NaN
index = find(isnan(tau_mean_fit)==1);
tau_mean_fit(index) = [];
tau_mean_broken_fit(index) = [];
index = find(isnan(tau_mean_broken_fit)==1);
tau_mean_fit(index) = [];
tau_mean_broken_fit(index) = [];

figure;
plot(tau_mean_broken_fit,tau_mean_fit,'.','MarkerSize',2);
p = polyfit(tau_mean_broken_fit,tau_mean_fit,1); 
f = polyval(p,tau_mean_broken_fit);
hold on;
plot(tau_mean_broken_fit,f,'-');

tauq_revised = (tau_mean_broken.*p(1))+p(2); 

figure(together);
hold all;
loglog(deltanq,tauq_revised.*1e6,'b--','LineWidth',3);

save('tau_surface.mat','deltanq','tauq_revised');

%% Given the surface lifetime, determine the SRV
load('tau_surface.mat'); 

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

