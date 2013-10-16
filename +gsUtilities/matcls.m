function results = matcls(Ainit, Y, P, R, verbo)
%
% Solve CLS with an unconstrained Newton method. Really it's an
% equality-constrained problem, but I've eliminated the constraints with an
% affine function.
%
% Ainit: initial network (n x n)
% Y: expression matrix (n x m)
% P: perturbation matrix (n x m)
% R: correction matrix (n x n)
% verbo: print optimization messages (default = false)
%
% results: structure with fields:
%   A: optimal network (n x n)
%   status: {'Solved', 'Failed'}
%

% Settings
alpha = 0.01;
beta = 0.5;
prec = 1e-6;
tol = eps;
max_steps = 100;
max_bt_steps = 100;
n = size(Ainit,1);
d = n*n;
verbose = false;
if (nargin > 4)
    verbose = verbo;
end
tab = char(9);

% If it's a zero network then don't waste your time
NonZeros = (abs(Ainit) >= tol);
k = sum(sum(NonZeros)); % dimension of the problem
if (k == 0)
    rstatus = 'Solved';
    results = struct;
    results.A = Ainit;
    results.status = rstatus;
    if (verbose)
        display(['Zero networks are trivial. Completed with status = "' rstatus '"']);
    end
    return;
end

% Build the affine matrix F
F = zeros(d,k);
nonzeros = reshape(NonZeros, d, 1);
ix = 1;
for i = 1:d
    if (nonzeros(i))
        F(i,ix) = 1;
        ix = ix + 1;
    end
end % initial point

% Presets
RRt = R*R';
YYt = Y*Y';
gterm2 = 2 .* RRt*P*Y';

% Check for rank-deficient Y
rankdefY = false;
m = size(Y, 2);
if (m < n)
    rankdefY = true;
end

% Hessian is constant
H = F' * kron(YYt, 2*RRt) * F;
if (rankdefY)
    Hinv = pinv(H); % pseudoinverse
else % Should be positive definite
    [H_L, notposdef] = chol(H);
    H_L = H_L'; % A = H_L * H_L'
    if (notposdef)
        error('Hessian is not positive definite');
    end
end

% Functions
function val = gfunc(X, Y, R, P)
    val = R'*X*Y + R'*P;
end
function val = f(X, Y, R, P)
    gx = gfunc(X, Y, R, P);
    val = trace(gx*gx');
end
function gradx = grad(A, RRt, YYt, gterm2, d)
    gradmatrix = 2 .* RRt*A*YYt + gterm2;
    gradx = reshape(gradmatrix, d, 1);
end

% Prepare
A = Ainit; % warm start
xhat = reshape(A, d, 1);
z = zeros(k,1); % the optimization variable
btx = 0;
dline = '--------------------------------------------------------------------------';
if (verbose)
    display(' ');
    display('OptiMatt CLS v2.0');
    display('------------------');
    str = sprintf('%u free variables, precision=%1.0e, alpha=%1.2f, beta=%1.2f',k,prec,alpha,beta);
    display([str ' ']);
    display([dline ' ']);
    display(['  step' tab '  value' tab tab '  gap' tab tab '  bt steps' tab '  || (AY+P)''*R ||']);
    display([dline ' ']);
end
rstatus = 'Failed';

% Newton's method for equality-constrained minimization
for i = 0:max_steps

    % Evaluate and get gradient
    fx = f(A, Y, R, P);
    g = F'*grad(A, RRt, YYt, gterm2, d);

    % Newton step and decrement
    if (rankdefY)
        znt = -Hinv*g;
    else
        znt = -H_L' \ (H_L \ g);
    end
    gap = -0.5*g'*znt;

    % Display progress
    if (verbose)
        frob = norm((A*Y+P)'*R, 'fro');
        if (i == 0)
            str = sprintf('  %u.\t  %+1.4e\t  %1.4e\t\t\t  %1.4e', i, fx, gap, frob);
            display([str ' ']);
        else
            str = sprintf('  %u.\t  %+1.4e\t  %1.4e\t  %u\t\t  %1.4e', i, fx, gap, btx, frob);
            display([str ' ']);
        end
    end

    % Check stopping criteria
    if (gap < prec)
        rstatus = 'Solved';
        break;
    end

    % Line search
    t = 1/beta;
    lhs = Inf;
    rhs = -Inf;
    btx = 0;
    while (true)
        t = beta*t;
        zcurr = z + t*znt;
        xcurr = F*zcurr + xhat;
        Acurr = reshape(xcurr, n, n);

        % Evaluate at new point to get lhs
        lhs = f(Acurr, Y, R, P);

        % Right-hand side
        rhs = fx + alpha*t*g'*znt;

        % Stop
        if (lhs < rhs)
            break;
        end

        % Short-circuit
        if (btx > max_bt_steps)
            error('Short-circuited line search');
        end
        btx = btx + 1;
    end

    % Update
    z = z + t*znt;
    x = F*z + xhat;
    A = reshape(x, n, n);

end % Newton steps

% Display result
if (verbose)
    display([dline ' ']);
    display(['Completed with status = "' rstatus '"']);
end

% Save results
results = struct;
results.A = A;
results.status = rstatus;

end