ndop = logspace(13,18,500); 

type = 'n'; 

deltan = ndop./1000; 

Nt = 1e12; 
Ect = 0;
Etv = 0.38;
sigma_n = 1.3e-14; 
sigma_p = 7e-17; 
T = 300;

for i = 1:length(ndop)
    [tau_SRH,p1,n1] = SRH_full_adv(Nt,sigma_n,sigma_p,Ect,Etv,T,deltan(i),ndop(i),type);
    tau_SRH_store(i) = tau_SRH; 
end

figure;
loglog(ndop,tau_SRH_store); 
%% SRH test for linearization - make "fake" data
clear all; close all; 
% dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\Method development\SRH test\Yan parameters';
deltan = logspace(12,18,500)'; 
Ndop = 1.9e16; 
resistivity = 0.529;
Nt = 1e12; 
Ect = 0.24;
Etv = 0;
sigma_n = 3.6e-14; 
sigma_p = 1e-14; 
T = 300;
type = 'p'; 
W = .02; 
SRV = 10.*ones(size(deltan)); %cm/s
tau_SRH_store = zeros(size(deltan)); 
tau_intr = zeros(size(deltan));
tau_surf = zeros(size(deltan)); 

for i = 1:length(deltan)
    [tau_SRH,p1,n1] = SRH_full_adv(Nt,sigma_n,sigma_p,Ect,Etv,T,deltan(i),Ndop,type);
    tau_SRH_store(i) = tau_SRH/1e6; 
    %Let's also get Richter lifetime
    tau_intr(i) = Richter(T,deltan(i),Ndop,type);
    D = diffusivity(T,type,Ndop,deltan(i));
    tau_surf(i) =(W/(2.*SRV(i)))+((1./D).*((W/pi)^2));
end
tau_meas = ((1./tau_SRH_store)+(1./tau_intr)+(1./tau_surf)).^(-1);

figure; 
loglog(deltan,tau_meas,'LineWidth',2); 
hold all;
loglog(deltan,tau_surf,'LineWidth',2);
hold all;
loglog(deltan,tau_intr,'LineWidth',2);
hold all;
loglog(deltan,tau_SRH_store,'LineWidth',2); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20); 
legend('Measured','Surface','Intrinsic','SRH'); 
title(['E_c-E_t = 0.23, k = ' num2str(sigma_n/sigma_p)],'FontSize',20);

lifetime_breakdown = struct('tau',tau_meas,'deltan',deltan,'tau_SRH',...
    tau_SRH_store,'tau_intr',tau_intr,'tau_surf',tau_surf);
save([dirname '\lifetime_breakdown_sim.mat'],'lifetime_breakdown');
SRVtoSave{1} = [deltan SRV];
save([dirname '\SRV_data.mat'],'SRVtoSave'); 

% dataSave{1} = [deltan,tau_meas];
% fileListShort{1} = 'Simulated_data_EcEt0-23';
% save([dirname '\Raw_data.mat'],'fileListShort','dataSave');
% info = struct('filename',fileListShort,'thickness',W,'resistivity',...
%     resistivity,'measured_resistivity',resistivity,'optical_constant',0.7,...
%     'calibration_factor',1,'temperature',T-273.15,'doping',Ndop);
% save([dirname '\meas_info.mat'],'info'); 
%% Make fake data, modifying input parameters. 
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\Method development\SRH test\modifying_raw';
deltan = logspace(12,18,500)'; 
Ndop = 1e16; 
rho = 1.46;
Nt = 1e12; 
Ect = 0;
Etv = 0.38;
sigma_n = 1.3e-14; 
sigma_p = 7e-17; 
T = 300;
type = 'p'; 
W = .02; 
Joe = 1e-15; 
tau_SRH_store = zeros(size(deltan)); 
tau_intr = zeros(size(deltan));
[tau_surf] = calculate_tau_surf_Joe(Joe,type,Ndop,deltan,W,T);
for i = 1:length(deltan)
    [tau_SRH,p1,n1] = SRH_full_adv(Nt,sigma_n,sigma_p,Ect,Etv,T,deltan(i),Ndop,type);
    tau_SRH_store(i) = tau_SRH/1e6; 
    %Let's also get Richter lifetime
    tau_intr(i) = Richter(T,deltan(i),Ndop,type);
end
tau_meas = ((1./tau_SRH_store)+(1./tau_intr)+(1./tau_surf)).^(-1);

figure; 
loglog(deltan,tau_meas,'LineWidth',2); 
hold all;
loglog(deltan,tau_surf,'LineWidth',2);
hold all;
loglog(deltan,tau_intr,'LineWidth',2);
hold all;
loglog(deltan,tau_SRH_store,'LineWidth',2); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20); 
legend('Measured','Surface','Intrinsic','SRH'); 
title(['E_t-E_v = 0.38, k = ' num2str(sigma_n/sigma_p)],'FontSize',20);

dataSave{1} = [deltan,tau_meas];
fileListShort{1} = 'Simulated_data_EtEv0-38_k185_raw';
thickness{1} = W; 
resistivity{1} = rho;
OC{1} = 0.7; 
calibration_factor{1} = 1; 
temperature{1} = T-273.15; 
doping{1} = Ndop; 

%Now we modify deltan + 10%
deltan_new = deltan.*1.1; 
dataSave{2} = [deltan_new,tau_meas];
fileListShort{2} = 'Simulated_data_EtEv0-38_k185_plusdeltan';
thickness{2} = W; 
resistivity{2} = rho;
OC{2} = 0.7; 
calibration_factor{2} = 1; 
temperature{2} = T-273.15; 
doping{2} = Ndop; 

%Modify deltan - 10%
deltan_new = deltan.*0.9; 
dataSave{3} = [deltan_new,tau_meas];
fileListShort{3} = 'Simulated_data_EtEv0-38_k185_minusdeltan';
thickness{3} = W; 
resistivity{3} = rho;
OC{3} = 0.7; 
calibration_factor{3} = 1; 
temperature{3} = T-273.15; 
doping{3} = Ndop; 

%Modify tau + 10%
tau_new = tau_meas.*1.1;  
dataSave{4} = [deltan,tau_new];
fileListShort{4} = 'Simulated_data_EtEv0-38_k185_plustau';
thickness{4} = W; 
resistivity{4} = rho;
OC{4} = 0.7; 
calibration_factor{4} = 1; 
temperature{4} = T-273.15; 
doping{4} = Ndop; 

%Modify tau -10%
tau_new = tau_meas.*0.9; 
dataSave{5} = [deltan,tau_new];
fileListShort{5} = 'Simulated_data_EtEv0-38_k185_minustau';
thickness{5} = W; 
resistivity{5} = rho;
OC{5} = 0.7; 
calibration_factor{5} = 1; 
temperature{5} = T-273.15; 
doping{5} = Ndop; 

save([dirname '\Raw_data.mat'],'fileListShort','dataSave');
info = struct('filename',fileListShort,'thickness',thickness,'resistivity',...
    resistivity,'measured_resistivity',resistivity,'optical_constant',OC,...
    'calibration_factor',calibration_factor,'temperature',temperature,'doping',doping);
save([dirname '\meas_info.mat'],'info'); 

%% PERC-specific
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\Method development\SRH test\modifying_raw\PERC\New files reattempt';
load([dirname '\meas_info.mat']); 
load([dirname '\Raw_data.mat']); 
saveName = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\Method development\SRH test\modifying_raw\PERC\New files reattempt\modifying param';
W = info(2).thickness; 
rho = info(2).resistivity; 
T = info(2).temperature+273.15; 
Ndop = info(2).doping; 
data = dataSave{2}; 
deltan = data(:,1); 
tau_meas = data(:,2); 

dataSave{1} = [deltan,tau_meas];
fileListShort{1} = 'Unaltered';
thickness{1} = W; 
resistivity{1} = rho;
OC{1} = 0.7; 
calibration_factor{1} = 1; 
temperature{1} = T-273.15; 
doping{1} = Ndop; 

%Now we modify deltan + 10%
deltan_new = deltan.*1.1; 
dataSave{2} = [deltan_new,tau_meas];
fileListShort{2} = 'Deltan plus ten';
thickness{2} = W; 
resistivity{2} = rho;
OC{2} = 0.7; 
calibration_factor{2} = 1; 
temperature{2} = T-273.15; 
doping{2} = Ndop; 

%Modify deltan - 10%
deltan_new = deltan.*0.9; 
dataSave{3} = [deltan_new,tau_meas];
fileListShort{3} = 'Deltan minus ten';
thickness{3} = W; 
resistivity{3} = rho;
OC{3} = 0.7; 
calibration_factor{3} = 1; 
temperature{3} = T-273.15; 
doping{3} = Ndop; 

%Modify tau + 10%
tau_new = tau_meas.*1.1;  
dataSave{4} = [deltan,tau_new];
fileListShort{4} = 'Tau plus ten';
thickness{4} = W; 
resistivity{4} = rho;
OC{4} = 0.7; 
calibration_factor{4} = 1; 
temperature{4} = T-273.15; 
doping{4} = Ndop; 

%Modify tau -10%
tau_new = tau_meas.*0.9; 
dataSave{5} = [deltan,tau_new];
fileListShort{5} = 'Tau minus ten';
thickness{5} = W; 
resistivity{5} = rho;
OC{5} = 0.7; 
calibration_factor{5} = 1; 
temperature{5} = T-273.15; 
doping{5} = Ndop; 

save([saveName '\Raw_data.mat'],'fileListShort','dataSave');
info = struct('filename',fileListShort,'thickness',thickness','resistivity',...
    resistivity','measured_resistivity',resistivity','optical_constant',OC',...
    'calibration_factor',calibration_factor','temperature',temperature','doping',doping');
save([saveName '\meas_info.mat'],'info'); 
%% SRH test for linearization - make "fake" data
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\Method development\SRH test';
deltan = logspace(12,18,500)'; 
Ndop = 1.9e16; 
resistivity = 0.529;
Nt = 1e12; 
Ect = 0.24;
Etv = 0;
sigma_n = 3.6e-14; 
sigma_p = 1e-14; 
T = [250 300 350 400];
type = 'p'; 
W = .02; 
SRV = 10.*ones(size(deltan)); %cm/s
tau_SRH_store = zeros(size(deltan)); 
tau_intr = zeros(size(deltan));
tau_surf = zeros(size(deltan)); 

for j = 1:length(T)
for i = 1:length(deltan)
    [tau_SRH,p1,n1] = SRH_full_adv(Nt,sigma_n,sigma_p,Ect,Etv,T(j),deltan(i),Ndop,type);
    tau_SRH_store(i) = tau_SRH/1e6; 
    %Let's also get Richter lifetime
    tau_intr(i) = Richter(T(j),deltan(i),Ndop,type);
    D = diffusivity(T(j),type,Ndop,deltan(i));
    tau_surf(i) =(W/(2.*SRV(i)))+((1./D).*((W/pi)^2));
end
tau_meas = ((1./tau_SRH_store)+(1./tau_intr)+(1./tau_surf)).^(-1);

figure; 
loglog(deltan,tau_meas,'LineWidth',2); 
hold all;
loglog(deltan,tau_surf,'LineWidth',2);
hold all;
loglog(deltan,tau_intr,'LineWidth',2);
hold all;
loglog(deltan,tau_SRH_store,'LineWidth',2); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20); 
legend('Measured','Surface','Intrinsic','SRH'); 
title(['E_c-E_t = 0.23, k = ' num2str(sigma_n/sigma_p) ', T = ' num2str(T(j))],'FontSize',20);

tau_meas_all{j} = tau_meas; 
tau_SRH_all{j} = tau_SRH_store; 
tau_surf_all{j} = tau_surf; 
tau_intr_all{j} = tau_intr; 
deltan_all{j} = deltan; 
SRVtoSave{j} = [deltan SRV];
dataSave{j} = [deltan,tau_meas];
fileListShort{j} = ['T' num2str(T(j)) 'C'];
thickness{j} = W; 
resistivity_store{j} = resistivity; 
optical_constant{j} = 0.7; 
doping{j} = Ndop; 
temp_store{j} = T(j)-273.15; 
end
lifetime_breakdown = struct('tau',tau_meas_all,'deltan',deltan_all,'tau_SRH',...
    tau_SRH_all,'tau_intr',tau_intr_all,'tau_surf',tau_surf_all);
save([dirname '\lifetime_breakdown_sim.mat'],'lifetime_breakdown');
save([dirname '\SRV_data.mat'],'SRVtoSave'); 

save([dirname '\Raw_data.mat'],'fileListShort','dataSave');
info = struct('filename',fileListShort,'thickness',thickness,'resistivity',...
    resistivity_store,'measured_resistivity',resistivity_store,'optical_constant',optical_constant,...
    'temperature',temp_store,'doping',doping);
save([dirname '\meas_info.mat'],'info'); 
