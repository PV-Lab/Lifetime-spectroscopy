function [tau_surf] = calculate_tau_surf_Joe(Joe,type,doping,deltan,W,T)
%Joe [A/cm2], type ['p' or 'n'], doping [cm^-3]',deltan vector to probe
%lifetime at [cm^-3], W sample thickness [cm],T sample temperature [K]

%Get the intrinsic carrier concentration for temperature
%Call standard density of states model
q = 1.602e-19; %C
%Boltzmann constant
k_B = 8.61733238e-5; %eV/K  
[NC,NV] = DOS_std(T); %cm^-3  
%Bandgap in silicon
if T == 300 || T == 298
    Eg = 1.1242; %eV
else
    [Eg] = Sze(T(tempi)); %eV
end
ni2 = NC*NV*exp(-Eg/(k_B*T)); %cm^-6

%Assume that Joe remains constant with temperature and calculate surface
%lifetime
tau_surf =(2.*Joe.*(doping+deltan)./(q*ni2*W)).^(-1); 