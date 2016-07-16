function MSE = calc_MSE(fit_parameters,X,tau)
%This function takes as input the fit parameters for a multi-defect fitting
%to actual lifetime and normalized injection level vectors. The output is to calculate
%the mean squared error of the fit. 
%fit_parameters should be a matrix where each row represents a different
%defect. Column 1 is the slope (m) of the linear fit, column 2 is the
%intercept (b). y = mx+b

[num_defects,num_parameters] = size(fit_parameters); 
%Check the size of num_parameters to be sure the input matrix is formatted
%correclty
if num_parameters~=2
    disp('There is an issue with the input matrix.'); 
end

tau_fit = zeros(size(tau)); 
for i = 1:num_defects
    %Calculate the lifetime associated with that defect
    tau_def = fit_parameters(i,1).*X+fit_parameters(i,2); %s
    %We actually need the inverse lifetime
    tau_def_contr = 1./tau_def; %s^-1
    %Add the contributino of the inverse lifetime to our total
    tau_fit = tau_fit+tau_def_contr; %s^-1
end
%We actually need to take the inverse of the fitted lifetime sum
tau_fit = 1./tau_fit; %s
%Now get the MSE
sum_sq = (abs(tau-tau_fit)).^2; 
MSE = sum(sum_sq)/length(sum_sq); 

    