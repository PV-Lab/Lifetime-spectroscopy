% Lifetime spectroscopy investigations
clear all; close all; 
%% Change in spectroscopy results with concentration - ptype
%First, what if we change the concentration (therefore the magnitude of the
%lifetime) of the same defect. k value ~ 26
directory = 'C:\Users\Mallory\Documents\Lifetime spectroscopy\IDLS investigations';
Ect = 0.4;
Etv = 0; 
Nt = [1e9 1e10 1e11 1e12];
sigma_n = 1.6e-14;
sigma_p = 6e-16;
type = 'p'; 
N_dop = 1e16; 
T = 300; 
deltan = logspace(12,18,500); 
deltan = deltan';
lifetime = figure; 
for i = 1:length(Nt)
    for j = 1:length(deltan)
        [tau_SRH,n1,p1] = SRH_full_std(Nt(i),sigma_n,sigma_p,Ect,Etv,T,deltan(j),N_dop,type);
        tau_SRH_store(j,i) = tau_SRH;
    end
    figure(lifetime); 
    loglog(deltan,tau_SRH_store(:,i),'-','LineWidth',2); 
    hold all; 
end
figure(lifetime); 
xlabel('Excess carrier density (cm^-^3)','FontSize',20);
ylabel('Lifetime (\mus)','FontSize',20);  
legend(num2str(Nt')); 
%Fit the "data"
fit_tries = 1e3; 
%Get sample parameters at specified temperature
[Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,N_dop,type); 
%This loop may be slow due to fitting
for i = 1:length(Nt)
    %Given the SRH lifetime, we just need to linearize the injection level
    %Plot the Murphy linearization for n-type
    %Normalized carrier density
    if type == 'p'
        X = (n0+deltan)./(p0+deltan);
    elseif type == 'n'
        X = (p0+deltan)./(n0+deltan);
    end
    X_store{i} = X; 
    [one_defect{i,1},MSE_one{i,1},two_defects{i,1},MSE_two{i,1},three_defects{i,1},MSE_three{i,1},all_parameters_store{i,1},all_MSE_store{i,1}] = fit_murphy_master(X,tau_SRH_store(:,i),T,directory,fit_tries);
end
%Now make the E-k curves for each set of "data"
defect_k = figure; 
defect_tau = figure; 
for i = 1:length(Nt)
    best_fit = one_defect{i,1}; 
    [m,n] = size(best_fit); 
    for j = 1:m
        [Et{i,j},k{i,j},alphanN{i,j}]=generate_Ek(best_fit(j,:),T,N_dop,type);
    end
    figure(defect_k); 
    plot(Et{i,1},k{i,1},'-','LineWidth',2); 
    hold all; 
    figure(defect_tau);
    plot(Et{i,1},1./alphanN{i,1},'LineWidth',2); 
    hold all; 
end
figure(defect_k); 
legend(num2str(Nt'));
%Plot true answer
hold all; 
plot([-1 2],[sigma_n/sigma_p sigma_n/sigma_p],'k--','LineWidth',2); 
hold all;
plot([1.124-Ect 1.124-Ect],[-100 100],'k--','LineWidth',2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20); 
figure(defect_tau); 
legend(num2str(Nt'));
hold all;
plot([1.124-Ect 1.124-Ect],[-100 100],'k--','LineWidth',2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
%% Change in spectroscopy results with concentration - ntype
%First, what if we change the concentration (therefore the magnitude of the
%lifetime) of the same defect. k value ~ 26
directory = 'C:\Users\Mallory\Documents\Lifetime spectroscopy\IDLS investigations';
Ect = 0.4;
Etv = 0; 
Nt = [1e9 1e10 1e11 1e12];
sigma_n = 1.6e-14;
sigma_p = 6e-16;
type = 'n'; 
N_dop = 1e16; 
T = 300; 
deltan = logspace(12,18,500); 
deltan = deltan';
lifetime = figure; 
for i = 1:length(Nt)
    for j = 1:length(deltan)
        [tau_SRH,n1,p1] = SRH_full_std(Nt(i),sigma_n,sigma_p,Ect,Etv,T,deltan(j),N_dop,type);
        tau_SRH_store(j,i) = tau_SRH;
    end
    figure(lifetime); 
    loglog(deltan,tau_SRH_store(:,i),'-','LineWidth',2); 
    hold all; 
end
figure(lifetime); 
xlabel('Excess carrier density (cm^-^3)','FontSize',20);
ylabel('Lifetime (\mus)','FontSize',20);  
legend(num2str(Nt')); 
%Fit the "data"
fit_tries = 1e3; 
%Get sample parameters at specified temperature
[Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,N_dop,type); 
%This loop may be slow due to fitting
for i = 1:length(Nt)
    %Given the SRH lifetime, we just need to linearize the injection level
    %Plot the Murphy linearization for n-type
    %Normalized carrier density
    if type == 'p'
        X = (n0+deltan)./(p0+deltan);
    elseif type == 'n'
        X = (p0+deltan)./(n0+deltan);
    end
    X_store{i} = X; 
    [one_defect{i,1},MSE_one{i,1},two_defects{i,1},MSE_two{i,1},three_defects{i,1},MSE_three{i,1},all_parameters_store{i,1},all_MSE_store{i,1}] = fit_murphy_master(X,tau_SRH_store(:,i),T,directory,fit_tries);
end
%Now make the E-k curves for each set of "data"
defect_k = figure; 
defect_tau = figure; 
for i = 1:length(Nt)
    best_fit = one_defect{i,1}; 
    [m,n] = size(best_fit); 
    for j = 1:m
        [Et{i,j},k{i,j},1./alphanN{i,j}]=generate_Ek(best_fit(j,:),T,N_dop,type);
    end
    figure(defect_k); 
    plot(Et{i,1},k{i,1},'-','LineWidth',2); 
    hold all; 
    figure(defect_tau);
    plot(Et{i,1},alphanN{i,1},'LineWidth',2); 
    hold all; 
end
figure(defect_k); 
legend(num2str(Nt'));
%Plot true answer
hold all; 
plot([-1 2],[sigma_n/sigma_p sigma_n/sigma_p],'k--','LineWidth',2); 
hold all;
plot([1.124-Ect 1.124-Ect],[-100 100],'k--','LineWidth',2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20); 
figure(defect_tau); 
legend(num2str(Nt'));
hold all;
plot([1.124-Ect 1.124-Ect],[-100 100],'k--','LineWidth',2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
%% Change in spectroscopy results with k value - ptype
directory = 'C:\Users\Mallory\Documents\Lifetime spectroscopy\IDLS investigations';
Ect = 0.4;
Etv = 0; 
Nt = [1e12];
sigma_n = [6e-17 6e-16 1.6e-14 1.2e-13];
sigma_p = 6e-16;
type = 'p'; 
N_dop = 1e16; 
T = 300; 
deltan = logspace(12,18,500); 
deltan = deltan';
lifetime = figure; 
for i = 1:length(sigma_n)
    for j = 1:length(deltan)
        [tau_SRH,n1,p1] = SRH_full_std(Nt,sigma_n(i),sigma_p,Ect,Etv,T,deltan(j),N_dop,type);
        tau_SRH_store(j,i) = tau_SRH;
    end
    figure(lifetime); 
    loglog(deltan,tau_SRH_store(:,i),'-','LineWidth',2); 
    hold all; 
end
figure(lifetime); 
xlabel('Excess carrier density (cm^-^3)','FontSize',20);
ylabel('Lifetime (\mus)','FontSize',20);  
legend(num2str((sigma_n./sigma_p)')); 
%Fit the "data"
fit_tries = 1e3; 
%Get sample parameters at specified temperature
[Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,N_dop,type); 
%This loop may be slow due to fitting
for i = 1:length(sigma_n)
    %Given the SRH lifetime, we just need to linearize the injection level
    %Plot the Murphy linearization for n-type
    %Normalized carrier density
    if type == 'p'
        X = (n0+deltan)./(p0+deltan);
    elseif type == 'n'
        X = (p0+deltan)./(n0+deltan);
    end
    X_store{i} = X; 
    [one_defect{i,1},MSE_one{i,1},two_defects{i,1},MSE_two{i,1},three_defects{i,1},MSE_three{i,1},all_parameters_store{i,1},all_MSE_store{i,1}] = fit_murphy_master(X,tau_SRH_store(:,i),T,directory,fit_tries);
end
%Now make the E-k curves for each set of "data"
defect_k = figure; 
defect_tau = figure; 
for i = 1:length(sigma_n)
    best_fit = one_defect{i,1}; 
    [m,n] = size(best_fit); 
    for j = 1:m
        [Et{i,j},k{i,j},alphanN{i,j}]=generate_Ek(best_fit(j,:),T,N_dop,type);
    end
    figure(defect_k); 
    plot(Et{i,1},k{i,1},'-','LineWidth',2); 
    hold all; 
    figure(defect_tau);
    plot(Et{i,1},1./alphanN{i,1},'LineWidth',2); 
    hold all; 
end
figure(defect_k); 
legend(num2str((sigma_n./sigma_p)'));
%Plot true answer
hold all;
plot([1.124-Ect 1.124-Ect],[-100 100],'k--','LineWidth',2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20); 
figure(defect_tau); 
legend(num2str((sigma_n./sigma_p)'));
hold all;
plot([1.124-Ect 1.124-Ect],[-100 100],'k--','LineWidth',2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
%% Change in spectroscopy results with k value - ntype
directory = 'C:\Users\Mallory\Documents\Lifetime spectroscopy\IDLS investigations';
Ect = 0.4;
Etv = 0; 
Nt = [1e12];
sigma_n = [6e-17 6e-16 1.6e-14 1.2e-13];
sigma_p = 6e-16;
type = 'n'; 
N_dop = 1e16; 
T = 300; 
deltan = logspace(12,18,500); 
deltan = deltan';
lifetime = figure; 
for i = 1:length(sigma_n)
    for j = 1:length(deltan)
        [tau_SRH,n1,p1] = SRH_full_std(Nt,sigma_n(i),sigma_p,Ect,Etv,T,deltan(j),N_dop,type);
        tau_SRH_store(j,i) = tau_SRH;
    end
    figure(lifetime); 
    loglog(deltan,tau_SRH_store(:,i),'-','LineWidth',2); 
    hold all; 
end
figure(lifetime); 
xlabel('Excess carrier density (cm^-^3)','FontSize',20);
ylabel('Lifetime (\mus)','FontSize',20);  
legend(num2str((sigma_n./sigma_p)')); 
%Fit the "data"
fit_tries = 1e3; 
%Get sample parameters at specified temperature
[Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,N_dop,type); 
%This loop may be slow due to fitting
for i = 1:length(sigma_n)
    %Given the SRH lifetime, we just need to linearize the injection level
    %Plot the Murphy linearization for n-type
    %Normalized carrier density
    if type == 'p'
        X = (n0+deltan)./(p0+deltan);
    elseif type == 'n'
        X = (p0+deltan)./(n0+deltan);
    end
    X_store{i} = X; 
    [one_defect{i,1},MSE_one{i,1},two_defects{i,1},MSE_two{i,1},three_defects{i,1},MSE_three{i,1},all_parameters_store{i,1},all_MSE_store{i,1}] = fit_murphy_master(X,tau_SRH_store(:,i),T,directory,fit_tries);
end
%Now make the E-k curves for each set of "data"
defect_k = figure; 
defect_tau = figure; 
for i = 1:length(sigma_n)
    best_fit = one_defect{i,1}; 
    [m,n] = size(best_fit); 
    for j = 1:m
        [Et{i,j},k{i,j},alphanN{i,j}]=generate_Ek(best_fit(j,:),T,N_dop,type);
    end
    figure(defect_k); 
    plot(Et{i,1},k{i,1},'-','LineWidth',2); 
    hold all; 
    figure(defect_tau);
    plot(Et{i,1},1./alphanN{i,1},'LineWidth',2); 
    hold all; 
end
figure(defect_k); 
legend(num2str((sigma_n./sigma_p)'));
%Plot true answer
hold all;
plot([1.124-Ect 1.124-Ect],[-100 100],'k--','LineWidth',2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20); 
figure(defect_tau); 
legend(num2str((sigma_n./sigma_p)'));
hold all;
plot([1.124-Ect 1.124-Ect],[-100 100],'k--','LineWidth',2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);

%% Change in SRH curves with temperature, different E values, asymmetric k
Ect = [0.1 0.4 0.7 0.9];
Etv = 0; 
Nt = 1e12;
sigma_n = 2.6667e-13;
sigma_p = 1e-14;
type = 'p'; 
N_dop = 1e16; 
T = [300 350 400 450]; 
deltan = logspace(12,18,500); 
deltan = deltan';
for i = 1:length(Ect)
    lifetime = figure; 
    for k = 1:length(T)
        for j = 1:length(deltan)
            [tau_SRH,n1,p1] = SRH_full_adv(Nt,sigma_n,sigma_p,Ect(i),Etv,T(k),deltan(j),N_dop,type);
            tau_SRH_store(j,k) = tau_SRH;
        end 
        loglog(deltan,tau_SRH_store(:,k),'-','LineWidth',2); 
        hold all; 
    end
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);  
    legend(num2str(T'));
    title(['E_c-E_t = ' num2str(Ect(i))]);
end
%% Change in SRH curves with temperature, different k values, constant E
Ect = [0.9];
Etv = 0; 
Nt = 1e12;
sigma_n = [6e-17 6e-16 1.6e-14 1.2e-13];
sigma_p = 6e-16;
k_value = sigma_n./sigma_p; 
type = 'p'; 
N_dop = 1e16; 
T = [300 350 400 450]; 
deltan = logspace(12,18,500); 
deltan = deltan';
for i = 1:length(k_value)
    lifetime = figure; 
    for k = 1:length(T)
        for j = 1:length(deltan)
            [tau_SRH,n1,p1] = SRH_full_adv(Nt,sigma_n(i),sigma_p,Ect,Etv,T(k),deltan(j),N_dop,type);
            tau_SRH_store(j,k) = tau_SRH;
        end 
        loglog(deltan,tau_SRH_store(:,k),'-','LineWidth',2); 
        hold all; 
    end
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);  
    legend(num2str(T'));
    title(['k = ' num2str(k_value(i))]);
end
%% Look at SRH curve along E-k curve
load('C:\Users\Mallory\Documents\Sinton visit\by sample\PERC May 26 revised resistivity\69-8\Ek_degraded_PVSC.mat'); 
Ect = 0; 
sigma_p = 1e-14; 
N_dop = 8.2e15; 
type = 'p'; 
Nt = 1e14; 
deltan = logspace(12,18,500); 
deltan = deltan'; 
indices = 1:1:length(Ek_degraded_PVSC); 
T = [25 100 175]; 
T = T+273.15; 
Et = Ek_degraded_PVSC(:,2); 
k_values = Ek_degraded_PVSC(:,1); 
for i = 1:length(indices); 
    lifetime = figure; 
    for k = 1:length(T)
        for j = 1:length(deltan)
            [tau_SRH,n1,p1] = SRH_full_adv(Nt,k_values(indices(i))*sigma_p,sigma_p,Ect,Et(indices(i)),T(k),deltan(j),N_dop,type);
            tau_SRH_store(j,k) = tau_SRH;
        end 
        loglog(deltan,tau_SRH_store(:,k),'-','LineWidth',2); 
        hold all; 
    end
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);  
    legend(num2str(T'));
    title(['k = ' num2str(k_values(indices(i))) ', Et = ' num2str(Et(indices(i)))]);
end


