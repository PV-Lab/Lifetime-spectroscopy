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

%%Load the data
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\2015 Oxygen-State Study\NOC samples\NOC 17 21 22 study\November 2015\One Day after Anneal\FZ';
process_xls_data(dirname);
%dataSave is a cell structure with probably one entry. The first column is
%injection level. The second column is lifetime. 

%% Given the measured lifetime, determine the SRV
% load('all_XLS_data.mat'); 
% data = dataSave{1}; %this assumes you're only looking at one sample
% deltan = data(:,1);
% tau = data(:,2);

load('C:\Users\Mallory\Documents\Australia\Passivation run\June 30 2016 passivation run\July 1 2016\16-6-28-P-2\averageTau.mat');
tau = tau_mean; 
deltan = deltanq;

N_dop = 5.7e15; %cm-3
W = 0.0280; %cm
T = 300; %K
type = 'p';

%Get the intrinsic lifetime
tau_intr = zeros(length(deltan),1);
%Calculate Richter model for the different doping levels
for i = 1:length(deltan)
    tau_intr(i,1) = Richter(T,deltan(i),N_dop,type);
end

%Get the surface lifetime
tau_surf = ((1./tau)-(1./tau_intr)).^(-1);
figure;
loglog(deltan,tau_surf);
hold all;
loglog(deltan,tau_intr);
hold all;
loglog(deltan,tau);
legend('Surface','Intrinsic','Measured');

D = 30.61; %minority carrier in n-type is hole 

SRV = W./((tau_surf-((1/D)*((W/pi)^2))).*2);

figure;
semilogx(deltan,SRV); 
xlabel('Excess carrier density [cm^{-3}]');
ylabel('SRV [cm/s]');

save('C:\Users\Mallory\Documents\Australia\Passivation run\June 30 2016 passivation run\July 1 2016\16-6-28-P-2\SRV_deltan.mat','deltan','SRV');

