function f = occupation_SRH(T,p,n,k_min,k_max,maxEc,minEv,Ec,Ev)
%Inputs:
%T in Kelvin
%p = total hole population
%n = total electron population
%k_min = minimum k value to make contour plot
%k_max = maximum k value to make contour plot
%Ec = conduction band location
%Ev = valence band location

%Get DOS at temperature based on advanced model
[NC, NV] = DOS_em(T);

%Boltzmann constant
k_B = 8.61733238e-5; %eV/K

%Bandgap in silicon, using Sze bandgap narrowing model
% [Eg] = Sze(T); %eV
% [Eg] = Sze(343);

%We are going to sweep through different k values and energy levels
k = linspace(k_min,k_max,10000); 
% Et = linspace(0,Eg,100); %eV, referenced to valence band edge
Et = linspace(minEv,maxEc,500);

f = zeros(length(k),length(Et)); 
[Eg_hold] = Sze(T); 

for i = 1:length(k)
    for j = 1:length(Et)
        if Et(j) > Ev && Et(j) < Ec
            Etv = Et(j)-Ev;
    %         Ect = Eg-Et(j);
    %         Ect = Eg_hold-Et(j);
            Ect = Ec-Et(j);
            %Calculate SRH densities
            n1 = NC*exp(-Ect/(k_B*T)); %cm^-3
            p1 = NV*exp(-Etv/(k_B*T)); %cm^-3
            f(i,j) = 1/(1+(((k(i)*n1)+p)/((k(i)*n)+p1))); 
    %         if Et(j)>Eg_hold
    %             f(i,j) = NaN;
    %         end
        else
            f(i,j) = NaN;
        end
    end
end

% figure;
% [X,Y] = meshgrid(Et',k); 
% surf(X,Y,f,'LineStyle','none');
% grid off;
% colormap('parula'); 
% view([0 90]); 
% xlabel('E_t-E_v [eV]','FontSize',20); 
% ylabel('k = \sigma_n/\sigma_p','FontSize',20); 
% axis([0 Eg k_min k_max]);

end