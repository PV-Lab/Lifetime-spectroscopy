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

function [deltan,tau] = remove_duplicates(deltan,tau)
%Deltan should be a vector. Tau can be a matrix but the number of rows
%should be equivalent to that of deltan. 
%We need to check and make sure no deltan (x-values) are repeated. This
%will interfere with the interpolation. 
    %Make the list where we will store indices to delete
    store_index = []; 
    %Initialize the count for the number of entries to be in our list
    count = 1; 
    %Loop through every entry and check for duplicates
    for j = 1:length(deltan)
        index = find(deltan==deltan(j)); 
        [m,n] = size(index); 
        %If there is more than one match, we need to delete an index. We
        %will just choose to delete the second index in the list. for normal
        %lifetime stuff this will be spurious low injection data anyway. 
        if m>1
            store_index(count:count+m-2) = index(2:end); 
            count = count+1; 
        end
    end
    deltan(store_index) = [];
    tau(store_index,:) = [];