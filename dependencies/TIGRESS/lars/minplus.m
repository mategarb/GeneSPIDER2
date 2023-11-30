function [m, I, J] = minplus(X)

% Find the minimum and its index over the (strictly) positive part of X matrix

% Remove complex elements and reset to Inf
[I,J] = find(0~=imag(X));
for i = 1:length(I),
    X(I(i),J(i)) = Inf;
end

X(X<=0) = Inf;
m = min(min(X));
[I,J] = find(X==m);

end