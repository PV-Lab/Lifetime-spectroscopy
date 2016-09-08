function [PC_corr] = correct_dark_voltage(time,PC,time_before)
time_before = time_before*0.75; 
%Now look for the points to find the dark voltage on PC
PC_dark = PC(find(time<=time_before)); 
PC_dark = nanmean(PC_dark); 
%Now correct all values for the dark voltage
PC_corr = PC-PC_dark; 
%Now replot the results
figure;
plot(time,PC);
hold all;
plot(time,PC_corr);
xlabel('Time [s]','FontSize',20); 
ylabel('Voltage [V]','FontSize',20);
legend('PC','PC corrected'); 
title('Dark voltage corrected','FontSize',20);