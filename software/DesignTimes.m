function fitted = DesignTimes(SVDAx, vecU)

% Computes A*vec(B) for the Gaussian solvers, where vecU holds the upper
% triangular entries of B and A is the design matrix after regressing out X
if strcmp(SVDAx.SP1_SolveMethod, 'svd')
    fitted = SVDAx.Vt'*(SVDAx.U'*vecU.*SVDAx.Sdiag);
else
    fitted = SVDAx.A'*vecU;
end

end
