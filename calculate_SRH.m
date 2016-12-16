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

function [tau_SRH,deltanq] = calculate_SRH(measured_file,SRV_file,deltanq,N_dop,T,type,W,D,saveFile,saveStart,sample_no)
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

% Interpolate the SRV so that it matches the measured lifetime
SRVq = interp1(deltan_SRV,SRV,deltanq);  

tau_surf =(W./(2.*SRVq))+((1/D).*((W/pi)^2)); %cm/s

for i = 1:length(deltanq)
    tau_intr(i,1) = Richter(T,deltanq(i),N_dop,type);
end

load(measured_file); 
for i = 1:length(dataSave)
    % deltan = data.deltanq;
    % tau = data.tau_mean; 
    % tau = tau';
    data = dataSave{i};
    deltan = data(:,1); 
    tau = data(:,2); 

    tauq = interp1(deltan,tau,deltanq);
%     tauq = tau; 

    h=figure;
    loglog(deltanq,tauq.*1e6);
    hold all; 
    loglog(deltanq,tau_surf.*1e6); 
    hold all; 
    loglog(deltanq,tau_intr.*1e6); 

    xlabel('Excess carrier density (cm^-^3)','FontSize',30);
    ylabel('Lifetime (\mus)','FontSize',30);
    axis([5e13 1e17 0 25000]);

    tau_SRH_hold = ((1./tauq)-(1./tau_intr)-(1./tau_surf)).^(-1);
    hold all;
    loglog(deltanq,tau_SRH_hold.*1e6);
    legend('Measured','Surface','Intrinsic','SRH');
    title(['Sample ' sample_no],'FontSize',30); 

    hgsave(h,[saveStart 'LifetimeBreakdown_' num2str(i)]);
    print(h,'-dpng','-r0',[saveStart 'LifetimeBreakdown_' num2str(i) '.png']); 
    tau_SRH(i,:) = tau_SRH_hold;
end
save(saveFile,'tau_SRH','tau_intr','tauq','tau_surf','deltanq');

