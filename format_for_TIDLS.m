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

%Read a .txt file which has been output from python and store the lifetime
%and carrier density. Put this into a format which will be accepted by
%TIDLS codes. 
function format_for_TIDLS(data_filename,saveDir,lifetime_type)
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
    if isempty(this_line_num) == 1
        %We either have the header or an empty line
        if line ~= 1
            %This is not the header, but signals that we are entering a new
            %line. Store the previous data and start a new matrix. 
            dataSave{sample_count} = data_now; 
            sample_count = sample_count+1; 
            %We only need to store the excess carrier density and lifetime
            new_line = 1; %signal that next time we make a new data_now vector
        end
    elseif isempty(this_line_num)==0
        if line == 2 || new_line == 1
            %We need to start the data_now storage
            if strcmp(lifetime_type,'PC')==1
                data_now = [this_line_num(1,2) this_line_num(1,4)];
            elseif strcmp(lifetime_type,'PL')==1
                data_now = [this_line_num(1,3) this_line_num(1,5)];
            else
                display('Issue with specifying lifetime measurement type.');
            end
            new_line = 0; 
        else
            %Then we just add this line to our existing structure
            if strcmp(lifetime_type,'PC')==1
                data_now(end+1,:) = [this_line_num(1,2) this_line_num(1,4)];
            elseif strcmp(lifetime_type,'PL')
                data_now(end+1,:) = [this_line_num(1,3) this_line_num(1,5)];
            else
                display('Issue with specifying lifetime measurement type.');
            end
        end
    end
    tline = fgetl(fileID);
    %Increment the line counter
    line = line+1; 
end
%Save the very last sample 
dataSave{sample_count} = data_now; 
%close the file
fclose(fileID); 
%Now we should have read the entire file. Save. 
save([saveDir '\Raw_data.mat'],'dataSave');  
