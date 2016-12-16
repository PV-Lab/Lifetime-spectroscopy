%{
MIT License

Copyright (c) [2016] [Mallory Ann Jensen, jensenma@alum.mit.edu]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%}

function [tau_mean,deltanq]=average_tau(filename,saveStart)
%average_tau(filename) processes a .mat file which has been produced by
%process_xls_data. All of the data sets within the file are assumed to
%belong the same sample. The curves are interpolated and then averaged to
%obtain an average curve.

load(filename); 

deltanq = logspace(10,17,1000); %typical max range for the QSSPC tool

together = figure; 

for i = 1:length(dataSave)
    data = dataSave{i};
    deltan = data(:,1); 
    tau = data(:,2);
    
    %remove spurious low injection data
    [deltan,tau] = remove_lowinj(deltan,tau,1e10);
    
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

print('-dpng','-r0',[saveStart '\AverageLifetime.png']);
saveFile = [saveStart '\averageTau.mat'];
save(saveFile,'tau_mean','deltanq');