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

function [tau_surf] = calculate_tau_surf_Joe(Joe,type,doping,deltan,W,T)
%Joe [A/cm2], type ['p' or 'n'], doping [cm^-3]',deltan vector to probe
%lifetime at [cm^-3], W sample thickness [cm],T sample temperature [K]

%Get the intrinsic carrier concentration for temperature
%Call standard density of states model
q = 1.602e-19; %C
%Boltzmann constant
k_B = 8.61733238e-5; %eV/K  
%Bandgap in silicon
if T <= 303 & T >= 297
    Eg = 1.1242; %eV
    %Also note that the background carrier concentration be assumed to be
    %equal to the doping
    N_dop = doping;
    [NC,NV] = DOS_std(T); %cm^-3  
else
    [Eg] = Sze(T); %eV
    [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,doping,type);
    if type == 'p'
        N_dop = p0;
    elseif type == 'n'
        N_dop = n0;
    end
    %Density of states effective mass model for consistency
    [NC,NV] = DOS_em(T); %cm^-3  
end
ni2 = NC*NV*exp(-Eg/(k_B*T)); %cm^-6

%Assume that Joe remains constant with temperature and calculate surface
%lifetime
tau_surf =(2.*Joe.*(N_dop+deltan)./(q*ni2*W)).^(-1); 