function A = scalefree(N,n,varargin)
% Create a scalefree network with N nodes and specific sparseness
% with preferential attachment
%
% A = datastruct.scalefree(N,n[, seed])
%
% N:    number of nodes
% n:    wished average number of links per node if n>=1, else relative sparsity of total possible links.
% pin:  The probability of adding an incoming edge to a node
% pout: pin+pout is the probabiltiy to add an outgoing edge to a given gene
%       pin and pout togheter controls the in/outgoing edges degree
%       in the powerlaw trying to ensure that the powerlaw is
%       primarily on outgoing edges. 
% seed: matrix to work as a seed, assumed to have size < N
%
% A:    undirected scalefree network matrix

% rng('shuffle');


rank_check = true;
pin = NaN;
pout = NaN;

if ~isempty(varargin)
    for i=1:length(varargin)
        if isa(varargin{i},'double')
            if isnan(pin)
                pin = varargin{i};
            elseif isnan(pout)
                pout = varargin{i};
            end
        elseif isa(varargin{i},'logical')
            rank_check = varargin{i};
        else
            seed = varargin{i};
        end
    end
end

% check that pin and pout is actually set as it will crash without this
if isnan(pin) 
    pin = 0.5;
end
if isnan(pout)
    pout = 0.3;
end

% set up the degree n, or relative sparsity S
if n < 1
    sparsity = n;
else
    sparsity = n/N;
end


m0 = round(sparsity*N); % recualculate back to a mean degree per node

if ~exist('seed','var')
    seed = logical(full(sprand(m0*2,m0*2-1,m0/m0^2)));
    k = 0;
    if rank_check
        while rank(double(seed)) < min(size(seed))
            seed = randn(size(seed)).*logical(full(sprand(m0*2,m0*2-1,m0/m0^2)));
            seed(floor(rand*numel(seed)+1)) = 1;
            k = k + 1;

            if ~mod(k,100)
                fprintf('k = %d\n',k)
            end
        end
    end
    tmp = zeros(m0*2);
    for i=1:size(seed,1)
        tmp(i,i+1:end) = seed(i,i:end);
    end

    for i=2:size(seed,1)
        tmp(i,1:i-1) = seed(i,1:i-1);
    end
    seed = logical(tmp); % creates small matrix that is initial for next steps
end

% check if there is any row with all zeros
% if there is, add 1 somewhere
zr = find(all(seed==0,2)); % rows
if ~isempty(zr)
    for z = 1:length(zr)
        seed(zr(z),randi([1,length(seed(zr(z),:))],[1,1])) = 1; % put inside a single 1 at random place
    end
end
% the same for columns
zc = find(all(seed==0,1)); % columns
if ~isempty(zc)
    for z = 1:length(zc)
        seed(randi([1,length(seed(zc(z),:))],[1,1]),zc(z)) = 1; % put inside a single 1 at random place
    end
end

% need to add checking if there are non-zero columns rows in seeds
A = zeros(N);
A(1:size(seed,1),1:size(seed,2)) = seed; % put initial matrix at the beginning of final matrix

% the part where the algorithm adds more links 
for i = (m0*2+1):N

    if rand < rem(sparsity, floor(sparsity)) % avoid rounding error
        m = ceil(m0);
    else
        m = floor(m0);
    end

    if m == 0
        m = 1;
    end

    k = 0;
    while k < m
        ps = 0;
        r = rand;
        for inode = 1:i-1

            % probability that a node will receive new edge
            pl = sum(abs(A(inode,1:i-1)))/nnz(A(1:i-1,1:i-1));

            ps = ps + pl; %nthroot(pl,N);% convert to biased probabilities, preserves single and double links to be in majority

            if r < ps % decide if the random probability is smaller than a given probability then add new edge
                
                r2 = rand();
                randC = rand;

                if randC <= 0.5 % randomize the sign of edge  50/50 chance, maybe it should be change cause it is more common to have positive regulation e.g. with TFs
                    val = -1;
                else
                    val = 1;
                end

                if r2 < pin % equal chance to be on upper and lower part, i.e. symetrical
                    A(i,inode) = val;
                elseif r2 < pin+pout
                    A(inode,i) = val;
                end
                
                k = k + 1;
                break
            end
        end
    end
end

% this part is added to reduce the extreme diagonalization that
% occurs when using +/-1 as the only values in the system. 
% A = A.*rand(N,N);
% this is added to ensure the stabilize script don't crash every
% other time it is used. 

% finnally we assume negative self loops on diagonal, otherwise NaNs occur after
% stabilization
A(eye(N)==1) = -1;

return