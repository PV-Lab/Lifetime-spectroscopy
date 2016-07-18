function [k] = calculate_Mod(T)
%This function calculates the k-value for a molybdenum donor defect at
%a given temperature (K). The fit is based on equations 6 and 7 in Paudyal JAP
%2010. Input can be a single temperature value or a vector of temperatures.
%The output structure will correspond to the input structure. 

for i = 1:length(T)
    k(i) = (6.9587e5)*(T(i)^-1.88); 
end

figure;
semilogy(T,k,'o','LineWidth',3); 
xlabel('Temperature [K]','FontSize',20); 
ylabel('k-value [-]','FontSize',20);
title('Mo_d','FontSize',20); 