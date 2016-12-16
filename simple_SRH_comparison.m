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

%Comparing SRH lifetimes
saveStart = 'C:\Users\Mallory\Documents\Sinton visit\May 1 2016\124-1';
filename1 = 'C:\Users\Mallory\Documents\Sinton visit\May 1 2016\124-1\lifetime_breakdown.mat'; 
filename2 = 'C:\Users\Mallory\Documents\Non-contact crucible\Australia experiments\June 22 2016\124-1\lifetime_breakdown.mat'; 

load(filename1); 

figure; 
loglog(deltanq,tau_SRH.*1e6);
tau_SRH1 = tau_SRH; 
load(filename2); 
hold all; 
loglog(deltanq,tau_SRH.*1e6);
tau_SRH2 = tau_SRH; 
legend('Measurement 1','Measurement 2'); 

xlabel('Excess carrier density (cm^-^3)','FontSize',30);
ylabel('Lifetime (\mus)','FontSize',30);
axis([5e13 1e17 0 10000]);
% print('-dpng','-r0',[saveStart '\SRH_comparison.png']);

figure;
semilogx(deltanq,((tau_SRH1-tau_SRH2)./tau_SRH1).*100); 