function [tau_mean,deltanq]=average_tau(filename)
%average_tau(filename) processes a .mat file which has been produced by
%process_xls_data. All of the data sets within the file are assumed to
%belong the same sample. The curves are interpolated and then averaged to
%obtain an average curve.

load(filename); 

deltanq = logspace(13,17,500); %typical max range for the QSSPC tool

together = figure; 

for i = 1:length(fileListShort)
    data = dataSave{i};
    deltan = data(:,1); 
    tau = data(:,2);
    
    %We need to check and make sure no deltan (x-values) are repeated. This
    %will interfere with the interpolation. 
    deltan = flipud(deltan);
    tau = flipud(tau); 
    [deltan,IX] = sort(deltan); 
    tau = tau(IX); 
    store_index = []; 
    count = 1; 
    for j = 1:length(deltan)
        index = find(deltan==deltan(j)); 
        [m,n] = size(index); 
        if m>1
            store_index(count:count+m-2) = index(2:end); 
            count = count+1; 
        end
    end
    deltan(store_index) = [];
    tau(store_index) = [];
    
    %Plot the resulting curve
    loglog(deltan,tau.*1e6,'.','MarkerSize',5);
    hold all;
    
    %Interpolate at the query points
    tauq(:,i) = interp1(deltan,tau,deltanq);
end

%Find the mean of the lifetime
tau_mean = nanmean(tauq,2); 

%Plot the result
loglog(deltanq,tau_mean.*1e6,'-','LineWidth',2); 

xlabel('Excess carrier density (cm^-^3)','FontSize',30);
ylabel('Lifetime (\mus)','FontSize',30);
set(gca,'FontSize',20);
set(gca,'LineWidth',2);

print('-dpng','-r0',['AverageLifetime.png']);