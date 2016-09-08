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