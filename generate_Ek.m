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

function [Et,k,alphanN]=generate_Ek(best_fit,T,doping,type)
%[Et,k,alphanN]=generate_Ek(best_fit,T,doping,type). This function takes as input a vector (one row, two columns) of the best
%fit for one defect. The first element is the slope of the fit, the second
%element is the intercept of the fit. T is temperature in Kelvin, doping is
%doping level, type is p or n. Both p-type and n-type analyses are coded.
%The fit parameters should be in the units of seconds for predictable
%units. If they are not, then the quantity alphanN will be scaled according
%to the input units. 

%Boltzmann constant
k_B = 8.61733238e-5; %eV/K  
[Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,doping,type); 
[NC,NV] = DOS_em(T); %cm^-3 
[Eg] = Sze(T); %eV
[vth_e,vth_h] = vth_em(T); %cm/s
%Define the energy levels for evaluation
%Let's try instead referencing to the intrinsic energy
Et_min = -Eiv; 
Et_max = Eg-Eiv; 
% Et = linspace(0,Eg,250); %eV
Et = linspace(Et_min,Et_max,250); %eV, Et-Ei
Q = zeros(size(Et)); 
alphanN = zeros(size(Et));
k = zeros(size(Et)); 
A = best_fit(1)+best_fit(2); %X -> 1
B = best_fit(2); %X -> 0
C = best_fit(1)/A; %slope/X -> 1
if type == 'p'
    for j = 1:length(Et)
        %Calculate n1
%         n1 = NC*exp(-(Eg-Et(j))/(k_B*T)); 
        n1 = NC*exp(-(Eg-Eiv-Et(j))/(k_B*T)); 
        %Calculate p1
%         p1 = NV*exp(-Et(j)/(k_B*T)); 
        p1 = NV*exp(-(Et(j)+Eiv)/(k_B*T)); 
        %Calculate the Q values for these defects
        Q(j) = (C+(p1/p0))/(1-(n1/p0)-C);
        %Calculate the quantity alphan*Nt for these defects
        alphanN(j) = (1/B)*(1+((1/p0)*((Q(j)*n1)+p1))); 
        %Calculate the k values for these defects
        k(j) = Q(j)*vth_h/vth_e; 
    end
elseif type == 'n'
    for j = 1:length(Et)
        %Calculate n1
%         n1 = NC*exp(-(Eg-Et(j))/(k_B*T)); 
        n1 = NC*exp(-(Eg-Eiv-Et(j))/(k_B*T)); 
        %Calculate p1
%         p1 = NV*exp(-Et(j)/(k_B*T)); 
        p1 = NV*exp(-(Et(j)+Eiv)/(k_B*T)); 
        %Calculate the Q values for these defects
        Q(j) = (1-(p1/n0)-C)/(C+(n1/n0)); 
        %Calculate the k values for these defects
        k(j) = Q(j)*vth_h/vth_e; 
        %Calculate the alphap*Nt for these defects
        alphanN(j) = (1/B)*(1+(n1/n0)+(p1/(n0*Q(j)))); 
    end
end

indices = find(k<0); 
k(indices) = NaN;
indices = find(alphanN<0); 
alphanN(indices) = NaN; 