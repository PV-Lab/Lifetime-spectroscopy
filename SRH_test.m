ndop = logspace(13,18,500); 

type = 'n'; 

deltan = ndop./1000; 

Nt = 1e12; 
Ect = 0;
Etv = 0.38;
sigma_n = 1.3e-14; 
sigma_p = 7e-17; 
T = 300;

for i = 1:length(ndop)
    [tau_SRH,p1,n1] = SRH_full_adv(Nt,sigma_n,sigma_p,Ect,Etv,T,deltan(i),ndop(i),type);
    tau_SRH_store(i) = tau_SRH; 
end

figure;
loglog(ndop,tau_SRH_store); 