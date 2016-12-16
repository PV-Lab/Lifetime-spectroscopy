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

function recommended_bias = bias_light_correction(filename,mode,deltan_input)
%This function takes as input a filename including the location of a Sinton
%WCT-120 file. Bias light corrections are applying to try to find an
%appropriate value which minimizes the effect of trapping. Mode is QSS,
%Generalized, or Transient. Deltan_input is the value that the lifetime
%should be evaluated at for the bias light. 

%We need to read in the minority carrier density as a function of the
%implied suns
data = xlsread(filename,'RawData','A4:E124'); 
time = data(:,1); 
ref_cell = data(:,3); 
deltan = data(:,4);
lifetime = data(:,5); 
%We also need the conversion factor V/sun
conversion = xlsread(filename,'Settings','C5'); 
%Convert measured voltage to suns
suns = ref_cell/conversion; %suns
bias_values = [0,round(min(suns)*100)/100:.01:3];
%We need to some sample parameters
data = xlsread(filename,'Summary','E2:K2');
thickness = data(1); 
optical_constant = data(3); 
des_lifetime = data(7); 

summary = figure;
loglog(deltan,lifetime,'k','LineWidth',2); 
for i = 1:length(bias_values)
    indices = find(suns>0.85*bias_values(i) & suns<1.15*bias_values(i)); 
    if isempty(indices) == 0
        %Now we need to plot the apparent carrier density versus suns
%         figure;
%         plot(suns(indices),deltan(indices),'.'); 
%         xlabel('Suns');
%         ylabel('Apparent carrier density [cm^-^3]'); 
%         title(['Bias light = ' num2str(bias_values(i))]);
        %We find the linear fit
        lin_fit = polyfit(suns(indices),deltan(indices),1); 
        %Evaluate our fit at the bias point
        deltan_bias = polyval(lin_fit,bias_values(i)); 
        %Subtract the apparent trap density from each value of deltan
        deltan_mod = deltan-deltan_bias; 
        %Calculate the equivalent traps
        traps = deltan_bias-(bias_values(i)*0.038*optical_constant*des_lifetime*1e-6/(1.6e-19*thickness)); 
        minority_carrier_density = deltan-traps; 
        %We also need to calculate the lifetime
    else 
        deltan_bias = 0; 
        minority_carrier_density = deltan; 
        deltan_mod = deltan; 
    end
    delta_suns = suns-bias_values(i); %suns
    generation = delta_suns.*0.038.*optical_constant./1.6e-19; %pairs/s
    if strcmp(mode,'QSS')==1
        lifetime_mod = deltan_mod.*thickness./generation; %s
    elseif strcmp(mode,'Generalized')==1
        for j = 1:length(deltan_mod)
            window_start = j-2; 
            if window_start <= 0
                window_start = 1; 
            end
            window_end = j+2; 
            if window_end > length(deltan_mod)
                window_end = length(deltan_mod); 
            end
            x = deltan_mod(window_start:window_end); 
            y = time(window_start:window_end); 
            fit = polyfit(x,y,1); 
            dndt(j,1) = fit(1); 
            if abs(dndt(j,1))<1000000000000
                dndt(j,1) = -1000000000000; 
            end
        end
        lifetime_mod = deltan_mod.*thickness./(generation-(thickness.*dndt)); 
    elseif strcmp(mode,'Transient')==1
        disp('No trapping correction can be made for transient mode.'); 
        lifetime_mod = deltan_mod.*thickness./(-thickness.*dndt); 
    end 
    %Now we plot the resulting lifetime
%     figure(summary); 
%     hold all;
%     loglog(minority_carrier_density,lifetime_mod,'--'); 
    %Grab the lifetime at our input deltan
    if deltan_input > deltan_bias
        %Get rid of any repeated values
        [minority_carrier_density,ia] = unique(minority_carrier_density,'stable');
        lifetime_mod = lifetime_mod(ia); 
        lifetime_store(1,i) = interp1(minority_carrier_density,lifetime_mod,deltan_input); 
    else
        disp('Bias input greater than injection level');
        lifetime_store(1,i) = NaN; 
    end
    lifetime_mod_store{i} = lifetime_mod; 
    deltan_mod_store{i} = minority_carrier_density; 
end
figure(summary); 
hold all; 
best_index = find(lifetime_store==min(lifetime_store)); 
loglog(deltan_mod_store{best_index},lifetime_mod_store{best_index},'--','LineWidth',2); 
legend('As measured','Corrected'); 
xlabel('Excess carrier density [cm^-^3]'); 
ylabel('Lifetime [s]'); 
figure; 
plot(bias_values,lifetime_store,'o'); 
xlabel('Bias light [suns]'); 
ylabel(['Measured lifetime at ' num2str(deltan_input) ' [s]']);
recommended_bias = bias_values(best_index); 
    
    
    