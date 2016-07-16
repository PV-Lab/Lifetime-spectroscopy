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
