%Average QSSPL data
function [dataSave] = average_QSSPL(dirname,sample,T,num_repeats,time_before,Ctonum)
    for i = 1:length(T)
        for j =1:num_repeats(i)
            %Make filename
%             filename = [dirname '\' sample '_' num2str(T(i)) Ctonum num2str(j)...
%                 '.Raw Data.dat'];
            filename = [dirname '\' sample Ctonum num2str(j)...
                '.Raw Data.dat'];
%             filename = [dirname '\' sample Ctonum num2str(j)...
%                 '.Raw Data.dat'];
%             filename = [dirname '\' sample '_' Ctonum num2str(j)...
%                 '.Raw Data.dat'];
            %Make the savename
            savename = [dirname '\' sample '_' num2str(T(i)) Ctonum num2str(j)...
                '_rawData.mat'];
            %Read the data
            [data_now,time,PC,PL,gen]=format_QSSPL_rawdata(filename,savename);
            %Correct the PC
            [PC_corr] = correct_dark_voltage(time,PC,time_before);
            %Let's assume the time vector is the same for all of the
            %measurements
            if j == 1
                PC_sum = PC_corr; 
                gen_sum = gen; 
                time_store = time; 
                PL_sum = PL; 
            else
                PC_sum = PC_sum+PC_corr; 
                gen_sum = gen_sum+gen; 
                PL_sum = PL_sum+PL; 
            end
            if isequal(time,time_store)==0
                disp('Error with time vectors, times are not equivalent sample-to-sample.'); 
            else
                time_store = time; 
            end
        end
        %Now we have the average for our given temperature
        PC_overall = PC_sum./num_repeats(i);
        gen_overall = gen_sum./num_repeats(i); 
        PL_overall = PL_sum./num_repeats(i); 
        %Plot the result
        figure;
        plot(time_store,gen_overall);
        hold all;
        plot(time,PL_overall); 
        hold all; 
        plot(time_store,PC_overall);
        xlabel('Time [s]','FontSize',20); 
        ylabel('Voltage [V]','FontSize',20);
        legend('Generation','PL','PC'); 
        title(['Average PC, T = ' num2str(T(i)) 'C'],'FontSize',20); 
        dataSave{i} = [time_store gen_overall PC_overall PL_overall]; 
    end