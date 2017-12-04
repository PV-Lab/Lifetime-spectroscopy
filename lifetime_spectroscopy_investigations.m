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
clear all; close all; 
directory = 'C:\Users\Mallory Jensen\Documents\Lifetime spectroscopy\IDLS investigations';
Ect = 0.4;
Etv = 0; 
Nt = [1e12];
sigma_n = [6e-20 6e-19 6e-18 6e-10];
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
Ect = [0.24];
Etv = 0; 
Nt = 1e12;
sigma_n = [3.6e-14];
sigma_p = 1e-14;
k_value = sigma_n./sigma_p; 
type = 'p'; 
N_dop = 1.9e16; 
T = [200 250 300 350 400]; 
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
% load('C:\Users\Mallory\Documents\Sinton visit\by sample\PERC May 26 revised resistivity\69-8\Ek_degraded_PVSC.mat'); 
Ect = 0; 
sigma_p = 1e-14; 
N_dop = 9.74e15; 
type = 'p'; 
Nt = 1e12; 
deltan = logspace(12,18,500); 
deltan = deltan'; 
% indices = 1:1:length(Ek_degraded_PVSC); 
T = [25 100 175]; 
T = T+273.15; 
% Et = Ek_degraded_PVSC(:,2); 
% k_values = Ek_degraded_PVSC(:,1); 
Et = [0.2 0.28 0.5 0.74 0.79];
k_values = [42.35 28.99 28.36 29.57 39.58];
indices = 1:1:length(Et); 
for i = 1:length(indices); 
    lifetime = figure; 
    for k = 1:length(T)
        for j = 1:length(deltan)
            [tau_SRH,n1,p1] = SRH_full_adv(Nt,k_values(indices(i))*sigma_p,sigma_p,Ect,Et(indices(i)),T(k),deltan(j),N_dop,type);
            tau_SRH_store(j,k) = tau_SRH;
        end 
        loglog(deltan,tau_SRH_store(:,k),'-','LineWidth',3); 
        hold all; 
    end
    xlabel('Excess carrier density (cm^-^3)','FontSize',30);
    ylabel('Lifetime (\mus)','FontSize',30);  
    legend(num2str((T-273.15)'));
    title(['k = ' num2str(k_values(indices(i))) ', E_t-E_v = ' num2str(Et(indices(i)))],'FontSize',30);
end
%% n-type variations with temperature
clear all; close all; clc;
directory = 'C:\Users\Mallory Jensen\Dropbox (MIT)\TIDLS data - UROP\Data\analysis by ingot\test trends'; 
type = 'n';
N_dop = 2e15; 
Nt = 1e12; 
T = [-75 -25 25]; T = T+273.15; 
Etv = [0.1 0.2 0.4 0.55 0.7 0.9 1.0];
k = [1e-4 1e-2 10]; 
sigma_p = 1e-15; 
sigma_n = k.*sigma_p; 
deltan = logspace(12,18,500); 
deltan = deltan';
Ect = 0;

for i = 1:length(Etv)
    for m = 1:length(sigma_n)
        lifetime = figure; 
        for k = 1:length(T)
            for j = 1:length(deltan)
                [tau_SRH,n1,p1] = SRH_full_adv(Nt,sigma_n(m),sigma_p,Ect,Etv(i),T(k),deltan(j),N_dop,type);
                tau_SRH_store(j,k) = tau_SRH;
            end 
            loglog(deltan,tau_SRH_store(:,k),'-','LineWidth',3); 
            hold all; 
        end
        xlabel('Excess carrier density (cm^-^3)','FontSize',30);
        ylabel('Lifetime (\mus)','FontSize',30);  
        legend(num2str((T-273.15)'));
        title(['k = ' num2str(sigma_n(m)/sigma_p) ', E_t-E_v = ' num2str(Etv(i))],'FontSize',30);
        tau_SRH_all{i,m} = tau_SRH_store; 
        h = gcf;
        hgsave(h,[directory '\tauSRH_k = ' num2str(sigma_n(m)/sigma_p) ', E_t-E_v = ' num2str(Etv(i)) '.fig']);
        print(h,'-dpng','-r0',[directory '\tauSRH_k = ' num2str(sigma_n(m)/sigma_p) ', E_t-E_v = ' num2str(Etv(i)) '.png']);
    end
end

%Fit the "data"
fit_tries = 1e3; 

for i = 1:length(Etv)
    for m = 1:length(sigma_n)
        tau_SRH_now = tau_SRH_all{i,m};
        defect_k = figure;
        defect_tau = figure;
        for k = 1:length(T)
            %Get sample parameters at specified temperature
            [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T(k),N_dop,type); 
            if type == 'p'
                X = (n0+deltan)./(p0+deltan);
            elseif type == 'n'
                X = (p0+deltan)./(n0+deltan);
            end 
            [one_defect,MSE_one,two_defects,MSE_two,three_defects,MSE_three,...
                all_parameters_store,all_MSE_store] ...
                = fit_murphy_master(X,tau_SRH_now(:,k),T(k),directory,fit_tries);

            %Now make the E-k curves for each set of "data"
            best_fit = one_defect;
            [Et,kval,alphanN]=generate_Ek(best_fit(1,:),T(k),N_dop,type);
            figure(defect_k); 
            semilogy(Et,kval,'-','LineWidth',2); 
            hold all; 
            figure(defect_tau);
            semilogy(Et,1./alphanN,'LineWidth',2); 
            hold all; 
        end
        figure(defect_k);
        xlabel('E_t-E_v [eV]','FontSize',20); 
        ylabel('k [-]','FontSize',20); 
        title(['k = ' num2str(sigma_n(m)/sigma_p) ', E_t-E_v = ' num2str(Etv(i))],'FontSize',30);
        legend(num2str(T'-273.15));
        hgsave(defect_k,[directory '\Ek_k = ' num2str(sigma_n(m)/sigma_p) ', E_t-E_v = ' num2str(Etv(i)) '.fig']);
        print(defect_k,'-dpng','-r0',[directory '\Ek_k = ' num2str(sigma_n(m)/sigma_p) ', E_t-E_v = ' num2str(Etv(i)) '.png']);
        figure(defect_tau);
        xlabel('E_t-E_v [eV]','FontSize',20); 
        ylabel('\tau_{n0} [s]','FontSize',20);
        title(['k = ' num2str(sigma_n(m)/sigma_p) ', E_t-E_v = ' num2str(Etv(i))],'FontSize',30);
        legend(num2str(T'-273.15));
        hgsave(defect_tau,[directory '\Etau_k = ' num2str(sigma_n(m)/sigma_p) ', E_t-E_v = ' num2str(Etv(i)) '.fig']);
        print(defect_tau,'-dpng','-r0',[directory '\Etau_k = ' num2str(sigma_n(m)/sigma_p) ', E_t-E_v = ' num2str(Etv(i)) '.png']);
    end
end