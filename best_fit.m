function [new_parameters_win,MSE_win,new_parameters_store,MSE_store]=best_fit(X,tau_SRH,start_values,fit_tries)

%Define the number of defects
[m,n] = size(start_values);
%Make a matrix of different fit parameters. Slope can be negative, but
%intercept cannot
start_values_matrix{1} = start_values; 
all_neg = start_values;
all_neg(:,1) = all_neg(:,1)*-1; 
start_values_matrix{2} = all_neg; 
count = 3;
for i = 1:m
    input = start_values; 
    input(i,1) = input(i,1)*-1;
    start_values_matrix{count} = input;
    count = count+1;
    for j = i+1:m
        input_now = input; 
        input_now(j,1) = input_now(j,1)*-1;
        start_values_matrix{count} = input_now; 
        count = count+1; 
    end
end

%Now try fminsearch for each of these combinations
for i = 1:length(start_values_matrix); 
    new_options = optimset('MaxFunEvals',fit_tries*m*n,'MaxIter',fit_tries*m*n); 
    [new_parameters_store{i},MSE_store(i),exitflag_store(i)] = fminsearch(@(fit_parameters) calc_MSE(fit_parameters,X,tau_SRH),start_values_matrix{i},new_options); 
end

%We want the minimum mean squared error
new_parameters_win = new_parameters_store{find(MSE_store==min(MSE_store))};
MSE_win = min(MSE_store); 
exitflag_win = exitflag_store(find(MSE_store==min(MSE_store))); 
if exitflag_win == 1
    disp_message = 'Converged';
elseif exitflag_win == 0
    disp_message = 'Maximum iterations/function evaluations reached.';   
elseif exitflag_win == -1
    disp_message = 'Algorithm terminated by output function.';
else
    disp_message = 'Unknown error.';
end
disp([num2str(m) ' defects - exit flag message: ' disp_message]);

%If there is a negative intercept on any of the fits, try this:
if ~isempty(find(new_parameters_win(:,2)<0)==1)
    disp('Trying a constrained fit because one of the y-intercepts was negative...'); 
    %We don't want a linear or nonlinear constraints
    A = []; b = []; Aeq = []; beq = [];
    [m,n] = size(start_values); 
    %We just want a lower bound = 0 on the intercept
    lb = [-Inf.*ones(m,1) zeros(m,1)]; 
    ub = Inf.*ones(m,n); 
    nonlcon=[];
%     myoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',fit_tries*m*n,'MaxIter',fit_tries*m*n);
    myoptions = optimoptions('fmincon','MaxFunEvals',fit_tries*m*n,'MaxIter',fit_tries*m*n);
    %Now try fmincon for each of these combinations
    original_length=(length(MSE_store));
    for i = 1:length(start_values_matrix); 
        [new_parameters,MSE,exitflag] = fmincon(@(fit_parameters) calc_MSE(fit_parameters,X,tau_SRH),start_values_matrix{i},A,b,Aeq,beq,lb,ub,nonlcon,myoptions);
        if MSE>0
            new_parameters_store{original_length+i}=new_parameters;
            MSE_store(original_length+i)=MSE;
            exitflag_store(original_length+i)=exitflag;
        end
    end
%     [new_parameters_store{length(new_parameters_store)+1},MSE_store(length(MSE_store)+1),exitflag_store(length(exitflag_store)+1)] = fmincon(@(fit_parameters) calc_MSE(fit_parameters,X,tau_SRH),start_values,A,b,Aeq,beq,lb,ub,nonlcon,myoptions);
    %We want the minimum mean squared error
    %Eliminate the negative intercept options
    search_fits = new_parameters_store;
    search_fits(find(MSE_store==MSE_win)) = []; 
    search_MSE = MSE_store;
    search_MSE(find(MSE_store==MSE_win))=[];
    new_parameters_win = search_fits{find(search_MSE==min(search_MSE))};
    MSE_win = min(search_MSE); 
    while ~isempty(find(new_parameters_win(:,2)<0)==1) 
        search_fits(find(search_MSE==MSE_win)) = [];
        search_MSE(find(search_MSE==MSE_win))=[];
        new_parameters_win = search_fits{find(search_MSE==min(search_MSE))};
        MSE_win = min(search_MSE); 
    end   
        
    exitflag_win = exitflag_store(find(MSE_store==min(MSE_store))); 
    if exitflag_win == 1
        disp_message = 'Converged';
    elseif exitflag_win == 0
        disp_message = 'Maximum iterations/function evaluations reached.';   
    elseif exitflag_win == -1
        disp_message = 'Algorithm terminated by output function.';
    else
        disp_message = 'Unknown error.';
    end
    disp(['Had to revise fit for negative intercept. ' num2str(m) ' defects - exit flag message: ' disp_message]);

end



