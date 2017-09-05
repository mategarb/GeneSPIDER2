function varargout = RSS(Alist,data,i,type)
% Function for calculating prediction Error PE for y and p,
%
% [RSSy RSSp] = RSS(Alist,data,i,type)
%
% If type='WRSS' is given then the weighted residual sum of squares is
% calculated.
% Otherwise the residual sum of squares is calculated.
%

if nargin < 4,
    type = 'RSS';
end

P = data.P;
Y = response(data);

for j = 1:size(Alist,3)
    A = Alist(:,:,j);
    if any(isnan(A))
        RSSy(j) = NaN;
        RSSp(j) = NaN;
        continue
    end
    y = -pinv(A)*P(:,i);
    p = -A*Y(:,i);

    if strcmp(type,'WRSS') % weighted residual sum of squares
        [cvY, cvP] = cov(data);
        if ~isempty(data.cvY), cvY = data.cvY; end
        if ~isempty(data.cvP), cvP = data.cvP; end
        if nnz(cvP) == 0, cvP = 1; end
        RSSy(j) = (Y(:,i) - y)'*pinv(cvY)*(Y(:,i) - y); % WRSS
        RSSp(j) = (P(:,i) - p)'*pinv(cvP)*(P(:,i) - p); % WRSS
    else
        RSSy(j) = norm(Y(:,i) - y, 'fro')^2;
        RSSp(j) = norm(P(:,i) - p, 'fro')^2;
    end
end

varargout{1} = RSSy;

if nargout > 1
    varargout{2} = RSSp;
end
