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

%Average QSSPL data
function [dataSave] = average_QSSPL(dirname,sample,T,num_repeats,time_before,Ctonum)
    for i = 1:length(T)
        for j =1:num_repeats(i)
            %Make filename
%             filename = [dirname '\' sample '_' num2str(T(i)) Ctonum num2str(j)...
%                 '.Raw Data.dat'];
%             filename = [dirname '\' sample Ctonum num2str(j)...
%                 '.Raw Data.dat'];
%             filename = [dirname '\' sample Ctonum num2str(j)...
%                 '.Raw Data.dat'];
            filename = [dirname '\' sample '_' Ctonum num2str(j)...
                '.Raw Data.dat'];
            %Make the savename
            savename = [dirname '\' sample '_' num2str(T(i)) 'C' Ctonum num2str(j)...
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