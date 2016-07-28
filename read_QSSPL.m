%Read and plot data from a QSSPL/PC measurement. The files should be .txt
%files with the header removed. 
function [deltanPL,deltanPC,tauPL,tauPC] = read_QSSPL(filename)
    data = textread(filename);
    time = data(:,1);
    deltanPL = data(:,2);
    deltanPC = data(:,2);
    deltanPL = data(:,3);
    tauPC = data(:,4);
    tauPL = data(:,5);
    figure;
    loglog(deltanPL,tauPL,'.');
    hold all;
    loglog(deltanPC,tauPC,'.');
    legend('PL','PC')
    xlabel('Excess carrier density [cm^-^3]','FontSize',20);
    ylabel('Lifetime [s]','FontSize',20); 
end