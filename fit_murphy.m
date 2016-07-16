function [p,order,MSE,MSE_store]=fit_murphy(X,tau_SRH,threshold)
%[p,order,MSE,MSE_store] = fit_murphy(X,tau_SRH,threshold): This function
%takes as input the normalized injection condition (X for p-type, Y for
%n-type) and the Shockley-Read-Hall lifetime. The lengths of the two
%vectors should be the same, and the entries should be correlated. This
%function attempts to fit linear lines to the data (one curve, room
%temperature) to derive information regarding defect parameters according
%to the publication by J. Murphy et al (Journal of Applied Physics 111,
%2012). One additional input is a threshold value for the mean squared
%error, used to identify a good fit to the data.

cpu_start =cputime; 

%Vary the initial fit parameters by... 
variation = 30; 

%Number of points to try fitting
N = 10; 

%First try a linear fit to the data
p = polyfit(X,tau_SRH,1); 
linear_fit = polyval(p,X); 

%Calculate mean-squared error
MSE = abs(tau_SRH-linear_fit); 
MSE = MSE.^2; 
MSE = (1/length(linear_fit))*sum(MSE); 

if MSE<=threshold
    order = 1; 
    %Plot the result
    figure; 
    h(1)=plot(X,tau_SRH); 
    hold all; 
    h(2)=plot(X,linear_fit); 
    legend(h,'Actual','Fit'); 
else
    %Start fit fitting the extremes to determine the range 
    X_new1 = X(1:5);
    tau_SRH_new1 = tau_SRH(1:5); 
    p_1 = polyfit(X_new1,tau_SRH_new1,1); 
    m1 = p_1(1); b1 = p_1(2); 
    if length(X)-14 > 0
        X_new2 = X(length(X)-14:end); 
        tau_SRH_new2 = tau_SRH(length(X)-14:end);
    else
        X_new2 = X(length(X)-2:end);
        tau_SRH_new2 = tau_SRH(length(X)-2:end);
    end
    p_2 = polyfit(X_new2,tau_SRH_new2,1);
    m2 = p_2(1); b2 = p_2(2); 
    
    slope_range_1 = m1-(m1*variation):2*variation*m1/N:m1+(m1*variation);
    slope_range_2 = m2-(m2*variation):2*variation*m2/N:m2+(m2*variation);
    intercept_range_1 = b1-(b1*variation):2*variation*b1/N:b1+(b1*variation);
    intercept_range_2 = b2-(b2*variation):2*variation*b2/N:b2+(b2*variation);
    
    
    MSE = 0; 
    count = 1; 
    %Try to fit two lines
    for i = 1:length(slope_range_1)
        for j = 1:length(intercept_range_1)
            for k = 1:length(slope_range_2)
                for m = 1:length(intercept_range_2)
                    
                    fit1 = (slope_range_1(i).*X)+intercept_range_1(j); 
                    fit2 = (slope_range_2(k).*X)+intercept_range_2(m); 
                    potential_fit = ((1./fit1)+(1./fit2)).^(-1); 
                    
                    MSE = abs(tau_SRH-potential_fit); 
                    MSE = MSE.^2; 
                    MSE = (1/length(linear_fit))*sum(MSE); 
                    
                    MSE_store(count,1)=slope_range_1(i);
                    MSE_store(count,2)=intercept_range_1(j); 
                    MSE_store(count,3)=slope_range_2(k); 
                    MSE_store(count,4)=intercept_range_2(m); 
                    MSE_store(count,5)=MSE;
                    
                    count = count+1; 
                end
            end
        end
    end
    
    MSE_values = MSE_store(:,5);
    index = find(MSE_values==min(MSE_values)); 
    MSE = MSE_store(index,5); 
    p = MSE_store(index,1:4); 
    order = 2; 
    
    %Plot the result
    fit1 = (p(1).*X)+p(2); 
    fit2 = (p(3).*X)+p(4); 
    potential_fit = ((1./fit1)+(1./fit2)).^(-1); 
%     figure; 
%     h(1)=plot(X,tau_SRH); 
%     hold all; 
%     h(2)=plot(X,fit1); 
%     hold all;
%     h(3)=plot(X,fit2);
%     hold all;
%     h(4)=plot(X,potential_fit); 
%     legend(h,'Actual','Fit 1','Fit 2','Harmonic Sum'); 
end

cpu_end = cputime; 
elapsed_time = cpu_end-cpu_start;

end
