function [one_defect,MSE_one,two_defects,MSE_two,three_defects,MSE_three,all_parameters_store,all_MSE_store] = fit_murphy_master(X,tau_SRH,T,directory,fit_tries)
%This function takes as input the normalized injection level and SRH
%lifetime and attempts to find the best fits for 1, 2, and 3 defects. 
    
    %First start with 1 defect
    p = polyfit(X,tau_SRH,1);
    [new_parameters_one,MSE_one,new_parameters_store_one,MSE_store_one]=best_fit(X,tau_SRH,p,fit_tries);
    %Find the value of the function at our new parameters
    linear_fit = polyval(new_parameters_one,X); 
    %Plot the result, saving space for the 1 and 2 defect results
    all_fits = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,3,1); 
    plot(X,tau_SRH,'b.'); 
    hold all; 
    plot(X,linear_fit,'k-','LineWidth',3); 
    legend('Measured','Fit'); 
    xlabel('X [-]','FontSize',15); 
    ylabel('Lifetime','FontSize',15); 
    %save the parameters and the MSE
    one_defect = new_parameters_one; 
    title({'One defect fit',['MSE = ' num2str(MSE_one)],['T = ' num2str(T) 'C']},'FontSize',25);
    
    %Now try 2 defects. First find starting parameters from fit_murphy function
    %Set the threshold at some arbitrary value
    threshold = 1e-15; 
    [p,order,MSE,MSE_store]=fit_murphy(X,tau_SRH,threshold); 
    if order == 2
        %Make the fit parameter matrix from the output p values
        fit_start = [p(1) p(2); p(3) p(4)];
    elseif order == 1
        %We need to pick some arbitrary start values
        X_new1 = X(1:5);
        tau_SRH_new1 = tau_SRH(1:5); 
        p_1 = polyfit(X_new1,tau_SRH_new1,1); 
        m1 = p_1(1); b1 = p_1(2); 
        X_new2 = X(length(X)-14:end); 
        tau_SRH_new2 = tau_SRH(length(X)-14:end);
        p_2 = polyfit(X_new2,tau_SRH_new2,1);
        m2 = p_2(1); b2 = p_2(2); 
        fit_start = [m1 b1; m2 b2];
    end
    %Find the best fit for two defects
    [new_parameters_two,MSE_two,new_parameters_store_two,MSE_store_two]=best_fit(X,tau_SRH,fit_start,fit_tries);
    %Plot the result
    figure(all_fits)
    subplot(1,3,2); 
    h(1)=plot(X,tau_SRH,'b.'); 
    hold all; 
    %Plot the result
    def1 = (new_parameters_two(1,1).*X)+new_parameters_two(1,2); 
    def2 = (new_parameters_two(2,1).*X)+new_parameters_two(2,2); 
    two_def_fit = ((1./def1)+(1./def2)).^(-1); 
    h(2)=plot(X,def1,'g--'); 
    hold all;
    h(3)=plot(X,def2,'r--');
    hold all;
    h(4)=plot(X,two_def_fit,'k-','LineWidth',3); 
    legend(h,'Actual','Fit 1','Fit 2','Harmonic Sum'); 
    xlabel('X [-]','FontSize',15); 
    ylabel('Lifetime','FontSize',15); 
    %save the parameters and the MSE
    two_defects = new_parameters_two; 
    title({'Two defects fit',['MSE = ' num2str(MSE_two)],['T = ' num2str(T) 'C']},'FontSize',25);
    
    %Finally, do this for three defects. Let's make some starting values
    fit_start_three = [one_defect;two_defects];
    %Find the best fit for two defects
    [new_parameters_three,MSE_three,new_parameters_store_three,MSE_store_three]=best_fit(X,tau_SRH,fit_start_three,fit_tries);
    %Plot the result
    figure(all_fits)
    subplot(1,3,3); 
    h(1)=plot(X,tau_SRH,'b.'); 
    hold all; 
    %Plot the result
    def1 = (new_parameters_three(1,1).*X)+new_parameters_three(1,2); 
    def2 = (new_parameters_three(2,1).*X)+new_parameters_three(2,2); 
    def3 = (new_parameters_three(3,1).*X)+new_parameters_three(3,2);
    three_def_fit = ((1./def1)+(1./def2)+(1./def3)).^(-1); 
    h(2)=plot(X,def1,'g--'); 
    hold all;
    h(3)=plot(X,def2,'r--');
    hold all;
    h(4)=plot(X,def3,'c--');
    hold all;
    h(5)=plot(X,three_def_fit,'k-','LineWidth',3); 
    legend(h,'Actual','Fit 1','Fit 2','Fit 3','Harmonic Sum'); 
    xlabel('X [-]','FontSize',15); 
    ylabel('Lifetime','FontSize',15); 
    %save the parameters and the MSE
    three_defects = new_parameters_three; 
    title({'Three defects fit',['MSE = ' num2str(MSE_three)],['T = ' num2str(T) 'C']},'FontSize',25);
    
    all_parameters_store = {new_parameters_store_one,new_parameters_store_two,new_parameters_store_three};
    all_MSE_store = {MSE_store_one,MSE_store_two,MSE_store_three}; 
    
    %We want to save the figure to look at later
    hgsave(all_fits,[directory '\Best fits_' num2str(round(T))]);
    print(all_fits,'-dpng','-r0',[directory '\Best fits_' num2str(round(T)) '.png']); 