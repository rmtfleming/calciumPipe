function [nb, elem] = find_nb_occurences(X)

% This function find the number of occurences of each element of a vector
% INPUT:
%     X: a vector
% OUTPUT:
%     nb: number of occurences of each element


nb = []; num = 100;

if length(X) > num   
    j = 1;
    for i = 1 : floor(length(X)/num)
        [n1, ~] = histc(X(i*num-num+1:i*num,1),...
            unique(X(i*num-num+1:i*num,1))); 
        nb = [nb; n1];
        % if the element nb 501 equals the previous one save the index to
        % merge later
        if X(i*num+1,1) == X(i*num,1)
            k(j) = length(nb);       % index of the elements to merge
            j=j+1;
        end
    end
    
    if length(X) - num*i > 0
        [n1, ~] = histc(X(num*i+1:end,1), unique(X(num*i+1:end,1))); 
        nb = [nb; n1];
    end    
end

% sum up the broken signals at nb 500
if exist('k', 'var')
    nb(k) = nb(k)+nb(k+1);
    nb(k+1) = [];
end

% vector of unique elements
for i = 1 : length(nb)
    elem(i) = X(sum(nb(1:i)));
end