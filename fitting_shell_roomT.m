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

clear all; close all; 
filename = 'C:\Users\Mallory\Documents\Australia\Quasi mono NTNU\QSSPC\A44\linSRH.mat';
directory = 'C:\Users\Mallory\Documents\Australia\Quasi mono NTNU\QSSPC\A44';
load(filename); 
[m,n] = size(tau_SRH); 
fit_tries = 1e5;
for i = 1:m
    tau_SRH_now = tau_SRH(i,:)'; 
    indices = find(tau_SRH_now<0 | isnan(tau_SRH_now)==1); 
    tau_SRH_now(indices) = [];
    X_now = X; 
    X_now(indices) = [];
    [one_defect,MSE_one,two_defects,MSE_two,three_defects,MSE_three,all_parameters_store,all_MSE_store] = fit_murphy_master(X_now,tau_SRH_now.*1e6,25,directory,fit_tries);
    to_write = zeros(6,3); 
    to_write(1:2,1) = one_defect';
    twodef = two_defects; 
    to_write(1:2,2) = twodef(1,:)';
    to_write(3:4,2) = twodef(2,:)';
    threedef = three_defects; 
    to_write(1:2,3) = threedef(1,:)';
    to_write(3:4,3) =threedef(2,:)';
    to_write(5:6,3) = threedef(3,:)';
    xlswrite([directory '\Linearized_data_wFit.xlsx'],[X_now,tau_SRH_now],['Sheet' num2str(i)]);
    xlswrite([directory '\Linearized_data_wFit.xlsx'],to_write,['Sheet' num2str(i)],'C1:E6'); 
end
%
%% Given best fit parameters, now create the Ek curves
clear all; close all; 
%Read in all the best fit parameters which have been sorted into similar
%defects
data = xlsread('C:\Users\Mallory\Documents\Australia\Quasi mono NTNU\QSSPC\fitting summary.xlsx','non grain boundary','B3:E10');
%indices in data that we want to plot together
matching_indices = [1 3;2 4;5 7;6 8]; 
doping = [1.5e16; 1.7e16; 1.6e16; 1.5e16];
type = 'p';
T = 300; 
doping_indices = [1 1 2 2 3 3 4 4];
[m,n] = size(matching_indices);
for i = 1:m
    def1_Ek = figure;
    def1_taun0 = figure;
    def2_Ek = figure;
    def2_taun0 = figure; 
    for j = 1:n
        %Make the plots for these pairs of indices
        doping_now = doping(doping_indices(matching_indices(i,j))); 
        for l = 1:2:3 %we have 2 defects
            best_fit = data(matching_indices(i,j),l:(l+1)); 
            [Et,k,alphaN] = generate_Ek(best_fit,T,doping_now,type);
            if l == 1
                figure(def1_Ek); 
                hold all; 
                plot(Et,k,'LineWidth',3); 
                figure(def1_taun0); 
                hold all; 
                plot(Et,1./alphaN,'LineWidth',3); 
            else
                figure(def2_Ek); 
                hold all;
                plot(Et,k,'LineWidth',3); 
                figure(def2_taun0); 
                hold all;
                plot(Et,1./alphaN,'LineWidth',3); 
            end
        end
    end
    figure(def1_Ek); 
    xlabel('E_t-E_v [eV]','FontSize',20); 
    ylabel('k [-]','FontSize',20);
    title(['Defect 1, Pair ' num2str(i)],'FontSize',20); 
    legend('Gettered','As-grown'); 
    figure(def1_taun0); 
    xlabel('E_t-E_v [eV]','FontSize',20); 
    ylabel('\tau_{n0} [s]','FontSize',20);
    title(['Defect 1, Pair ' num2str(i)],'FontSize',20); 
    legend('Gettered','As-grown');
    figure(def2_Ek); 
    xlabel('E_t-E_v [eV]','FontSize',20); 
    ylabel('k [-]','FontSize',20);
    title(['Defect 2, Pair ' num2str(i)],'FontSize',20); 
    legend('Gettered','As-grown');
    figure(def2_taun0); 
    xlabel('E_t-E_v [eV]','FontSize',20); 
    ylabel('\tau_{n0} [s]','FontSize',20);
    title(['Defect 2, Pair ' num2str(i)],'FontSize',20); 
    legend('Gettered','As-grown');
end
    
