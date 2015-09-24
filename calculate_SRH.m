function [tau_SRH,deltanq] = calculate_SRH(measured_file,SRV_file,deltanq,N_dop,T,type,W,D)
%[tau_SRH,deltan1] =
%calculate_SRH(measured_file,SRV_file,deltanq,N_dop,T,type,W,D). Measured
%file and SRV file should be of the format: 'filename.mat.' deltanq is
%logspace injection range targeted. N_dop is the doping level of the sample
%of interest in cm^-3. T is the temperature of the measured in K. type is
%either 'n' or 'p'. W is the thickness of the sample in cm. D is the
%diffusivity of the minority carrier in this sample in cm/s which can be
%obtained from PVCDROM. 


load(SRV_file);
deltan_SRV = deltan;

%Interpolate the SRV so that it matches the measured lifetime
SRVq = interp1(deltan_SRV,SRV,deltanq); 

tau_surf =(W./(2.*SRVq))+((1/D).*((W/pi)^2)); %cm/s

for i = 1:length(deltanq)
    tau_intr(i,1) = Richter(T,deltanq(i),N_dop,type);
end

load(measured_file); 
[m,n] = size(dataSave); 

for i = 1:m

    data = dataSave{i}; 
    deltan = data(:,1); 
    tau = data(:,2); 

    %Interpolate the lifetime data to use with the SRV later
    tauq = interp1(deltan,tau,deltanq); 

    figure;
    loglog(deltanq,tauq.*1e6);
    hold all; 
    loglog(deltanq,tau_surf.*1e6); 
    hold all; 
    loglog(deltanq,tau_intr.*1e6); 

    xlabel('Excess carrier density (cm^-^3)','FontSize',30);
    ylabel('Lifetime (\mus)','FontSize',30);
    axis([5e13 1e17 0 25000]);

    tau_SRH_hold = ((1./tau_meas)-(1./tau_intr)-(1./tau_surf)).^(-1);
    hold all;
    loglog(deltanq,tau_SRH_hold.*1e6);
    legend('Measured','Surface','Intrinsic','SRH');
    title(['Sample ' num2str(i)],'FontSize',30); 
    
    tau_SRH(:,i) = tau_SRH_hold; 
    
end