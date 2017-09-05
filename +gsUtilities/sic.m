function varargout = sic(Y,A)
% Calculates the strong irrepresentability condition of Y
% Zhao and Yu 2006 (http://portal.acm.org/citation.cfm?id=1248547.1248637)
%
% [irr, SIC] = sic(Y,A)
%
% Y: NxM expression matrix, N genes, M samples.
% A: network matrix.
%
% irr = min(1-SIC(~logical(A)))
% SIC = abs(Phiz'*Phipc*pinv(Phipc'*Phipc)*sign(A(i,A(i,:)~=0))'), if A(i,j) = 0
%

[N,M] = size(A);
Phi = Y';
nPhiz = [];
nPhipc = [];
reg0 = [];
regN0 = [];
for i = 1:N
    Phiz = Phi(:,A(i,:) == 0);
    Phipc = Phi(:,A(i,:) ~= 0);

    nPhiz = [nPhiz, Phiz(:)'];
    nPhipc = [nPhipc, Phipc(:)'];

    reg0tmp = [];
    for k=1:size(Phiz,2)
        reg0tmp(k) = norm(Phiz(:,k));
    end
    reg0 = [reg0, reg0tmp];

    regN0tmp = [];
    for k=1:size(Phipc,2)
        regN0tmp(k) = norm(Phipc(:,k));
    end
    regN0 = [regN0, regN0tmp];

    lregDiff(i) = min(regN0tmp)/max(reg0tmp);

    sic = abs(Phiz'*Phipc*pinv(Phipc'*Phipc)*sign(A(i,A(i,:)~=0))');
    k = 1;
    for j=1:N
        if A(i,j) == 0
            SIC(i,j) = sic(k);
            k = k + 1;
        end
    end
end

% SIC = SIC;
irr = min(1-SIC(~logical(A)));
if nargout == 0
    % irr = sum(sum(SIC < 1))/data.N^2;
    varargout{1} = irr;
end

if nargout >= 1
    % irr = sum(sum(SIC < 1))/data.N^2;
    varargout{1} = irr;
end

if nargout >= 2
    varargout{2} = SIC;
end

if nargout >= 3
    varargout{3} = nPhiz;
end

if nargout >= 4
    varargout{4} = nPhipc;
end

if nargout >= 5
    varargout{5} = reg0;
end

if nargout >= 6
    varargout{6} = regN0;
end

if nargout >= 7
    varargout{7} = lregDiff;
end
