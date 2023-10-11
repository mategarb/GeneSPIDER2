function [unique_degree, freq_degree, mean_degree] = genes_degree(A, io)
%GENES_DEGREE Summary of this function goes here
% A - network

% reading from rows to columns, as outdegree
if isequal(io,'out')
    dgr = sum(abs(logical(A)),2);
elseif isequal(io,'in')
    dgr = sum(abs(logical(A)),1);
elseif isequal(io,'all')
    dgr = sum(abs(logical(A)),1) + sum(abs(logical(A)),2)';
end
mean_degree = mean(dgr);

unique_degree = unique(dgr);

freq_degree = [];
for i = 1:length(unique_degree)
    freq_degree(i) = length(find(dgr == unique_degree(i)));
end

end

