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

%Read a .txt file which has been output from a measurement and store the
%RAW DATA. Put this into a format which will be accepted by other scripts.
function [data_now,time,PC,PL,gen] = format_QSSPL_rawdata(data_filename,savename)
%Open the text file
fileID = fopen(data_filename); 
%Read line-by-line 
tline = fgetl(fileID);
line = 1; 
sample_count = 1; 
new_line = 0; 
while ischar(tline)
    %If we have numbers, str2num should give us a 1D array
    this_line_num = str2num(tline); 
    %Check whether this is the header. If not, proceed. 
    if isempty(this_line_num)==0
        if line == 2 || new_line == 1
            %We need to start the data_now storage
            data_now = [this_line_num(1,1) this_line_num(1,2) this_line_num(1,3) this_line_num(1,4)];
            new_line = 0; 
        else
            %Then we just add this line to our existing structure
            data_now(end+1,:) = [this_line_num(1,1) this_line_num(1,2) this_line_num(1,3) this_line_num(1,4)];
        end
    end
    tline = fgetl(fileID);
    %Increment the line counter
    line = line+1; 
end
%close the file
fclose(fileID); 
%Plot our results for the visual check
figure;
plot(data_now(:,1),data_now(:,2));
hold all;
plot(data_now(:,1),data_now(:,3));
hold all;
plot(data_now(:,1),data_now(:,4));
xlabel('Time [s]','FontSize',20); 
ylabel('Voltage [V]','FontSize',20);
legend('Generation','PC','PL'); 
time = data_now(:,1); 
PC = data_now(:,3); 
gen = data_now(:,2);
PL = data_now(:,4); 
dataSave{1} = data_now; 
%Now we should have read the entire file. Save. 
save(savename,'dataSave','PC','time','PL','gen');  