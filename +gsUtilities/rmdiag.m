function varargout = rmdiag(net)
% removes diagonal elements and shifts the matrix along the second dimension.
%
% A = rmdiag(A)
%

if isa(varargin{i},'GeneSpider.Network') || isa(varargin{i},'tools.NetworkComparison')
    A = net.A;
elseif isa(varargin{i},'double')
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
