function [B_est, beta_est, count] = ProxyLogistic(SVDAx, X, y, Dk, Wk, delta, tol, max_iter, B0, beta0)

% IRLS solver for matrix-on-scalar logistic regression with Frobenius regularization% Inputs:
%   US      - n x n the product of U and S taken from SVD decompositon of A
%   V1      - n x p^2 right eigenvectors of A (only spinning the kernel space of A)
%   X       - n x d matrix of demographic covariates
%   y       - n x 1 binary response vector
%   delta   - regularization parameter
%   tol     - convergence tolerance
%   max_iter - maximum number of IRLS iterations 
%   v1       = V1'*vecY , for Y being a given matrix
%   v2       = (eye(p^2) - V1*V1')*vecY , for Y being a given matrix
%   B0      - (optional) if provided, then B0 is used for initialization of B
%   beta0   - (optional) if provided, then beta0 is used for initialization of beta
% SVDAx.SP1_SolveMethod decides how each IRLS step is solved:
%   'svd'      - US and V1 are used, n x n system
%   'woodbury' - SVDAx.A (n x pp design matrix) and SVDAx.K = A*A' are used,
%                n x n system
%   'direct'   - SVDAx.A is used, pp x pp system

% Outputs:
%   B_est   - estimated p x p matrix B
%   beta_est - estimated d x 1 vector beta
%   count   - number of IRLS loops
% -------------------------------------------------------------------------

% create Y 
Y = Dk + Wk/delta;

if SVDAx.UseSymmetricityAndZeroDiag
    delta = 2*delta;
end

% denote SVD matrices
p = SVDAx.pOriginal;
n = size(X,1);
useSVD = strcmp(SVDAx.SP1_SolveMethod, 'svd');
useWoodbury = strcmp(SVDAx.SP1_SolveMethod, 'woodbury');
if useSVD
    V1 = SVDAx.V;
    US = SVDAx.US;
    pp = size(V1,1); % pp = p*(p-1)/2 if UseSymmetricityAndZeroDiag == true, pp = p^2 if UseSymmetricityAndZeroDiag == false
else
    A = SVDAx.A;
    pp = size(A,2);
end

if nargin > 8
    % vecB0 = reshape(B0, [], 1);
    vecB = reshape(B0, [], 1);
    if SVDAx.UseSymmetricityAndZeroDiag
        vecB = vecB(SVDAx.idxs);
    end
    if useSVD
        AvecB = US*(V1'*vecB);
    else
        AvecB = A*vecB;
    end
else
    vecB = zeros(pp, 1);
    AvecB = zeros(size(X,1),1);
end

if nargin > 9
    beta = beta0;
else
    beta = zeros(size(X,2), 1);
end

% create two auxilliary vectors
vecY = reshape(Y, [], 1);
diagY = diag(Y);
if SVDAx.UseSymmetricityAndZeroDiag
    vecY = vecY(SVDAx.idxs);
end
if useSVD
    v1 = SVDAx.V'*vecY;
    V1t_vecY = V1'*vecY;
    V1V1t_vecY = V1*V1t_vecY;
    v2 = vecY - V1V1t_vecY;
    % v2 = (eye(pp) - V1*V1')*vecY;
elseif useWoodbury
    K = SVDAx.K;
    AvecY = A*vecY;
end

count = 0;
for iter = 1:max_iter
    eta = AvecB + X*beta;
    mu = 1 ./ (1 + exp(-eta));
    w = mu .* (1 - mu);
    z = eta - (mu - y) ./ (w + eps); % adjusted response
    
    % update z_tilde
    if useSVD
        WX = w.*X; % w.*X = W*X
        W2 = diag(w) - WX*((X'*WX)\WX');
        Q1 = US'*W2*US + delta*eye(n);
        z_tilde = Q1\(US'*W2*z + delta*v1);

        % Generates new estimates vectB and beta
        vecB_new = V1*z_tilde + v2;
        AvecB = US*z_tilde;
        beta_new = (X'*WX)\(WX'*(z-AvecB));
    elseif useWoodbury
        % beta and vecB = vecY + A'*alpha, for alpha, come from the n x n
        % system with R = delta*I + S*K*S, S = diag(sqrt(w))
        s = sqrt(w);
        R = (s*s').*K;
        R(1:n+1:end) = R(1:n+1:end) + delta;
        Xs = s.*X;
        rs = s.*(z - AvecY);
        [L, cholFlag] = chol(R, 'lower');
        if cholFlag == 0
            Rinv_rs = L'\(L\rs);
            Rinv_Xs = L'\(L\Xs);
        else
            % R is positive definite unless the weights are not finite
            Rinv_rs = R\rs;
            Rinv_Xs = R\Xs;
        end
        beta_new = (Xs'*Rinv_Xs)\(Xs'*Rinv_rs);
        alpha = s.*(Rinv_rs - Rinv_Xs*beta_new);
        vecB_new = vecY + A'*alpha;
        AvecB = AvecY + K*alpha;
    else
        % Solves (A'*W2*A + delta*I)*vecB = A'*W2*z + delta*vecY directly
        WX = w.*X; % w.*X = W*X
        W2 = diag(w) - WX*((X'*WX)\WX');
        Q1 = A'*W2*A + delta*eye(pp);
        vecB_new = Q1\(A'*W2*z + delta*vecY);
        AvecB = A*vecB_new;
        beta_new = (X'*WX)\(WX'*(z-AvecB));
    end

    % Check convergence
    converged = norm(vecB_new - vecB) + norm(beta_new - beta) < tol;
    vecB = vecB_new;
    beta = beta_new;
    if converged
        break;
    else
        count = count + 1;
    end
end

% Reshape and collect B and beta

% Reconstruct the symmetric matrix
if SVDAx.UseSymmetricityAndZeroDiag
    B_est         = zeros(p^2, 1);
    B_est(SVDAx.idxs,:) = vecB; 
    B_est         = reshape(B_est, [p, p]); % p-by-p matrix
    B_est         = B_est + B_est';
    B_est         = B_est + diag(diagY);
else
    B_est = reshape(vecB, p, p);
end
beta_est = beta;

end