function [possible_MCD,Joe]=Joe_curve(filename,min_MCD,max_MCD,doping)
%This function takes as input the filename (+location) of a new Sinton
%QSSPC spreadsheet and the minimum and maximum carrier densities. The
%function obtains a curve of Joe by manipulating the Excel spreadsheet -
%writing MCD values and then reading the result. The last input is the
%doping level which is used for plotting a reference value for the Joe
%curves. 
%Define the range of possible MCDs
possible_MCD = linspace(min_MCD,max_MCD,40);
for i = 1:length(possible_MCD)
    %Write the MCD to the file
    xlswrite(filename,possible_MCD(i),'User','F6');
    %Now that we've written it, go in a read the Joe
    Joe(i) = xlsread(filename,'Summary','L2'); 
end
%Plot the results
figure;
plot(possible_MCD,Joe,'o');
hold all; 
plot([doping doping],[0 max(Joe)],'k-');
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('J_o_e [A/cm^2]','FontSize',20); 
