function [k] = calculate_Tidd(T)
%This function calculates the k-value for a titanium double donor defect at
%a given temperature (celcius). The fit is based on Figure 7 in Paudyal JAP
%2009. Input can be a single temperature value or a vector of temperatures.
%The output structure will correspond to the input structure. 

for i = 1:length(T)
    k(i) = 0.39*T(i)+29.88; 
end

figure;
semilogy(T,k,'o','LineWidth',3); 
xlabel('Temperature [K]','FontSize',20); 
ylabel('k-value [-]','FontSize',20);
title('Ti_d_d','FontSize',20); 
