% Proximity operator for function F, where
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%             F(B) := sum_i (y_i - <A_i, B>)^2                        %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A - assumed to be 3-way tensor
% y - vector of observations
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------
%         Author:    Damian Brzyski
%         Date:      April 24, 2018
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set C instead of D for pure nuclear norm

% it solves problem 


%  argmin_B {  sum_i (y_i - <A_i, B>)^2 + delta || B - Dk - Wk/delta ||_F^2  }

function Bnew = ProxFsvd(y, SVDAx, Dk, Wk, delta)

%% 
p = size(Wk,1);

%%
idxs   = SVDAx.idxs;

%%
Ddelta     = Dk + Wk/delta;
DdeltaVec  = reshape(Ddelta, [p^2,1]);
DdeltaVecU = DdeltaVec(idxs);
ydelta     = y - DesignTimes(SVDAx, DdeltaVecU);

%% Ridge solution
persistent idChol deltaChol L % Cholesky factor, computed again only for a new design matrix or delta
switch SVDAx.SP1_SolveMethod
    case 'svd'
        AVty       = SVDAx.Vt * ydelta;
        MAVty      = (AVty.*SVDAx.Sdiag)./(SVDAx.Sdiag.^2 + 2*delta); % element-wise
        BridgeVec  = SVDAx.U * MAVty;
    case 'woodbury'
        if isempty(idChol) || idChol ~= SVDAx.id || deltaChol ~= delta
            R         = SVDAx.K;
            R(1:size(R,1)+1:end) = R(1:size(R,1)+1:end) + 2*delta;
            L         = chol(R, 'lower');
            idChol    = SVDAx.id;
            deltaChol = delta;
        end
        BridgeVec  = SVDAx.A*(L'\(L\ydelta));
    case 'direct'
        if isempty(idChol) || idChol ~= SVDAx.id || deltaChol ~= delta
            R         = SVDAx.AAt;
            R(1:size(R,1)+1:end) = R(1:size(R,1)+1:end) + 2*delta;
            L         = chol(R, 'lower');
            idChol    = SVDAx.id;
            deltaChol = delta;
        end
        BridgeVec  = L'\(L\(SVDAx.A*ydelta));
end

%% Reconstruct the symmetric matrix
Bridge         = zeros(p^2, 1);
Bridge(idxs,:) = BridgeVec; 
Bridge         = reshape(Bridge, [p, p]); % p-by-p matrix
Bridge         = Bridge + Bridge';

%% B update
Bnew           = Bridge + Ddelta;

end
