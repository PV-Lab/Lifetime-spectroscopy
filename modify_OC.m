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

%Modify the optical constant directly in the files. 
function [deltan_trans,tau_trans,dataSave_gen,OC_values,MSE,min_index]=modify_OC(filename_trans,filename_gen,OC_range,inj_range,resolution)
    %Read the transient data which is our reference
    [deltan,tau] = read_lifetime_data(filename_trans);
    %Save the raw data
    saveRaw_trans = [deltan,tau]; 
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,inj_range(2));
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,inj_range(1));
    deltan_trans = deltan_rev; tau_trans = tau_rev; 
    %Cycle through the possible optical constants
    OC_values = linspace(OC_range(1),OC_range(2),resolution); 
    for possible_value = 1:length(OC_values)
        %First we write the possible value
        xlswrite(filename_gen,OC_values(possible_value),'User','E6');
        %Now we read what we have written
        [deltan,tau] = read_lifetime_data(filename_gen);
        saveRaw_gen{possible_value} = [deltan,tau]; 
        [deltan_rev,tau_rev] = remove_highinj(deltan,tau,inj_range(2));
        [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,inj_range(1));
        %Now interpolate at the reference injection values
        [tau_rev] = interp1(deltan_rev,tau_rev,deltan_trans); 
        %Save the data
        dataSave_gen{possible_value} = [deltan_trans,tau_rev];
        %Let's grab the MSE
        %Now get the MSE
        sum_sq = (abs(tau_rev-tau_trans)).^2; 
        nan_index = find(isnan(sum_sq)==1); 
        if isempty(nan_index)==0
            sum_sq(nan_index)=[]; 
        end
        MSE(possible_value) = sum(sum_sq)/length(sum_sq); 
    end
    %Now plot the MSE versus the optical constant
    figure;
    plot(OC_values,MSE,'o','MarkerSize',8,'LineWidth',2); 
    xlabel('Optical constant [-]','FontSize',20); 
    ylabel('Mean squared error','FontSize',20); 
    %Find the minimum mean squared error
    min_index = find(MSE == min(MSE)); 
    %Plot the lifetime result for the minimum MSE
    figure;
    loglog(saveRaw_trans(:,1),saveRaw_trans(:,2)); 
    hold all;
    dataSave_win = saveRaw_gen{min_index};
    loglog(dataSave_win(:,1),dataSave_win(:,2)); 
    xlabel('Excess carrier density [cm^-^3]','FontSize',20); 
    ylabel('Lifetime [s]','FontSize',20); 
    title('Fit with minimum mean squared error','FontSize',20); 
    disp(['The winning optical constant is...' num2str(OC_values(min_index))]); 
end
    
        
