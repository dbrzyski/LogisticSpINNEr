% Proximity operator for function G, where
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                  G(C) := lambda_N*||C||_*                            %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ||C||_* - nuclear norm of C
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------
%         Author:    Damian Brzyski
%         Date:      April 24, 2018
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
function Cnew = ProxG(D, W2, delta, lambda_N, solOptions)

% solOptions (optional) - singleThreadCalc: if 1, eig runs on one thread
%                       - eigSymmetric: if true, symmetric form of eig and
%                         reconstruction only from the kept eigenvalues
if nargin < 5
    solOptions = struct('singleThreadCalc', 0, 'eigSymmetric', false);
end

Ddelta    = D + W2/delta;
Ddelta    = (Ddelta + Ddelta')/2;
if solOptions.singleThreadCalc == 1
    nt = maxNumCompThreads(1);
end
if solOptions.eigSymmetric
    [U, S]    = eig(Ddelta);  % Ddelta = U*S*U';
else
    [U, S, V] = eig(Ddelta);  % Ddelta = U*S*V';
end
if solOptions.singleThreadCalc == 1
    maxNumCompThreads(nt);
end
diagS     = diag(S);
if solOptions.eigSymmetric
    sTsh      = sign(diagS).*max(abs(diagS) - lambda_N/delta, 0); % soft tresholoding of singular values
    kept      = sTsh ~= 0;
    Cnew      = U(:, kept)*(sTsh(kept).*U(:, kept)');
else
    Stsh      = diag(  sign(diagS).*max(abs(diagS) - lambda_N/delta, 0)  ); % soft tresholoding of singular values
    Cnew      = U*Stsh*V';
end

end
