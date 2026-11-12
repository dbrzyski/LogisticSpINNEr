function out = minNormEstim(y, SVDAx)

%% Objects
p = SVDAx.pOriginal;

idxs   = SVDAx.idxs;

%% Estimate
switch SVDAx.SP1_SolveMethod
    case 'svd'
        % Regressing out X leaves singular values equal to zero up to rounding
        nonzero = SVDAx.Sdiag > max(size(SVDAx.U))*eps(max(SVDAx.Sdiag));
        Vty    = SVDAx.Vt(nonzero,:) * y./SVDAx.Sdiag(nonzero);
        Bvec   = SVDAx.U(:,nonzero)*Vty;
    case 'woodbury'
        Bvec   = SVDAx.A*(pinv(SVDAx.K)*y);
    case 'direct'
        Bvec   = pinv(SVDAx.A')*y;
end
Bestim         = zeros(p^2, 1);
Bestim(idxs,:) = Bvec; 
Bestim         = reshape(Bestim, [p, p]); % p-by-p matrix
Bestim         = Bestim + Bestim';
out            = struct;
out.B          = Bestim;

end
