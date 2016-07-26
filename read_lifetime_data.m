%[deltan,tau] = read_lifetime_data(filename). Given a single filename, this
%function reads the lifetime data. 
function [deltan,lifetime] = read_lifetime_data(filename)
    %For the NEW spreadsheet format
    data = xlsread(filename,'RawData','E4:G124');
    lifetime = data(:,1); %seconds
    deltan = data(:,3); %cm^-3
    figure;
    loglog(deltan,lifetime,'.');
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (seconds)','FontSize',20);
end