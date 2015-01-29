function varargout = rmdiag(net)
% Removes diagonal elements and shifts upper triangular elements -1 along the second dimension.
% Ad = rmdiag(A)
%
% A: initial network (n x n)
% Ad: non diag network (n x n-1)
%

if isa(net,'GeneSpider.Network') || isa(net,'tools.NetworkComparison')
    A = net.A;
elseif isa(net,'double')
    A = net;
end

tmp = [];
for i=1:length(A(1,1,:))
    Atmp = A(:,:,i);
    l = tril(Atmp,-1);
    l(:,end) = [];
    u = triu(Atmp,1);
    u(:,1) = [];
    tmp(:,:,i) = u+l;
end

varargout{1} = tmp;
