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

function [DPSSk, DPSStaun0, Ect] = DPSS(T,N_dop,type,deltan_fit,tau_fit)
%[DPSSk,DPSStaun0,Ect] = DPSS(T,N_dop,type,deltan_fit,tau_fit): This
%function takes as input one temperature value, the doping level of the
%sample, the doping type, and one set of experimental injection/lifetime
%measurements (SRH should already be separate). The output of this function
%is the solution surface for k and tau_n0 as a function of energy for this
%set of data. This function should be used in a loop to produce the set of
%DPSS for TIDLS.

%Define the fermi level and the equilibrium carrier concentrations
[Efi,Efv,p0,n0] = adv_Model_gen(T,N_dop,type); 

%Bandgap in silicon, using Sze bandgap narrowing model
[Eg] = Sze(T); %eV

%Define the temperature-dependent effective densities of state
[NC, NV] = DOS_em(T);

kB = 8.61733238e-5; %eV/K

%Define the fit type.
tau = fittype( 'taun0*(((p0+(NV*exp((Ect-Eg)/(kB*T)))+deltan)/(p0+n0+deltan))+(k*((n0+(NC*exp(-Ect/(kB*T)))+deltan)/(p0+n0+deltan))))','dependent',{'tau'},'independent',{'deltan'},'problem',{'p0','n0','kB','Eg','Ect','T','NC','NV'});
options = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0 0]);
Ect = 0:Eg/100:Eg; %Make solutions for defects placed throughout the bandgap
DPSSk = zeros(size(Ect));
DPSStaun0 = zeros(size(Ect));

for i = 1:length(Ect)
    
    %Do the fit and get the results
    result = fit(deltan_fit,tau_fit,tau,options,'problem',{p0,n0,kB,Eg,Ect(i),T,NC,NV});
    c = coeffvalues(result);
    DPSSk(i) = c(1);
    DPSStaun0(i) = c(2); 
 
end

end