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
%This script does a set of simple functions to perform lifetime fitting for
%measurements taken throughout degradation

function [easy_summary,all_defect] = fit_procedure(figure_handle,deltan,...
    tau,savename,T,N_dop,type,cutoff_low,cutoff_high)
    figure(figure_handle); 
    %First ask the user how to crop the data in high then low injection
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff_high);
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff_low);
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    %Save this figure for future reference if needed
    hgsave(figure_handle,savename);
    print(figure_handle,'-dpng','-r0',[savename '.png']);
    %Get the relevant carrier densities so that we can normalize
    [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,N_dop,type); 
     %Normalized carrier density
    if type == 'p'
        X = (n0+deltan_rev)./(p0+deltan_rev);
    elseif type == 'n'
        X = (p0+deltan_rev)./(n0+deltan_rev);
    end
    [two_defects,MSE_two,all_parameters_store,all_MSE_store] = fit_murphy_two(X,tau_rev,T,savename,1e5);
    [m,n] = size(two_defects);
    Et = cell(size(two_defects)); k = cell(size(two_defects)); 
    alphanN = cell(size(two_defects)); 
    for x = 1:m
        [Et{x},k{x},alphanN{x}]=generate_Ek(two_defects(x,:),T,N_dop,type);
    end
    all_defect = {two_defects,MSE_two,Et,k,alphanN,X,tau_rev}; 
    %Pick the dominant defect in low injection (3e14 to be safe)
    inj = find(abs(deltan_rev-3e14)==min(abs(deltan_rev-3e14))); 
    actual = tau_rev(inj); 
    def1 = two_defects(1,1)*X(inj)+two_defects(2,2);
    def2 = two_defects(2,1)*X(inj)+two_defects(2,2);
    if abs(def1-actual)<abs(def2-actual)
        Et_now = Et{1}; 
        Et_index = find(abs(0-Et_now)==min(abs(0-Et_now))); 
        k_now = k{1}; 
        alphanN_now = alphanN{1}; 
        %Defect 1 is dominant
        easy_summary = [k_now(Et_index) alphanN_now(Et_index)];
        %Defect 2 is secondary
        Et_now = Et{2}; 
        Et_index = find(abs(0-Et_now)==min(abs(0-Et_now))); 
        k_now = k{2}; 
        alphanN_now = alphanN{2};
        easy_summary(3)=k_now(Et_index);
        easy_summary(4)=alphanN_now(Et_index);
        easy_summary(5)=MSE_two; 
    elseif abs(def2-actual)<abs(def1-actual)
        Et_now = Et{2}; 
        Et_index = find(abs(0-Et_now)==min(abs(0-Et_now))); 
        k_now = k{2}; 
        alphanN_now = alphanN{2}; 
        %Defect 2 is dominant
        easy_summary = [k_now(Et_index) alphanN_now(Et_index)];
        %Defect 1 is secondary
        Et_now = Et{1}; 
        Et_index = find(abs(0-Et_now)==min(abs(0-Et_now))); 
        k_now = k{1}; 
        alphanN_now = alphanN{1};
        easy_summary(3)=k_now(Et_index);
        easy_summary(4)=alphanN_now(Et_index);
        easy_summary(5)=MSE_two; 
    end
end