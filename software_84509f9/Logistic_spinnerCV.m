%-------------------------------------------------------------------------- 
% This function finds the optimal pair of regularization parameters, 
% lambda_N and lambda_L, and finds SpINNEr estimates by solving the
% problem:
%--------------------------------------------------------------------------
%    argmin_{B, beta} {   0.5*sum_i ( y_i - X*beta - <A_i, B> )^2 + 
%                                 lambda_N*|| B ||_* + 
%                               lambda_L*|| vec(B o W) ||_1    }
%--------------------------------------------------------------------------
% In order to do so, two-dimensional cross-validation is applied and the 
% specific implementation of ADMM algorithm is used

% The SpINNEr toolbox is free software: you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% The SpINNEr toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------
%         Authors:    Damian Brzyski and Xixi Hu
%         Date:       June 27, 2018
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%-------------------------------------------
% REQUIRED INPUTS ARGUMENTS:               -
%-------------------------------------------
% y:                 The n-dimensional vector of responses
%~~~~~~~~~~~~~~~~~~~~~
% X:                 matrix of covariates included in the model but not
%                    penalized. Symbol "[]" should be provided to omit it.
%~~~~~~~~~~~~~~~~~~~~~
% AA:                This should be three-dimensional array with dimensions
%                    p,p and n, where ith slice, i.e. A(:,:,i), is a 
%                    symmetric matrix with zeros on the diagonal.
%                    The alternative form, AA = [A1, A2, A3,...,An],is also
%                    supported
%
%-------------------------------------------
% OPTIONAL INPUT ARGUMENTS (provided in standard name-value pair fashion)
%-------------------------------------------
% 'UseParallel':    Turns on/off the parallel computation for 
%                   cross-validation. Possible falues are: false (default)
%                   and true
%~~~~~~~~~~~~~~~~~~~~~
% 'W':              Matrix of weights in the SpINNEr optimization problem.
%                   Symmetric p-by-p matrices with nonnegative entries are
%                   possible. The default W has zeros on the diagonal and
%                   ones on the off-diagonal entries
%~~~~~~~~~~~~~~~~~~~~~
%'gridLengthN'      The number of different values of lambda_N considered
%                   in the cross-validation. These values are automatically
%                   selected.
%~~~~~~~~~~~~~~~~~~~~~
%'gridLengthL'      The number of different values of lambda_L considered
%                   in the cross-validation. These values are automatically
%                   selected.   
%~~~~~~~~~~~~~~~~~~~~~
%'kfolds'           The number of folds considered in cross-validation
%~~~~~~~~~~~~~~~~~~~~~
%'displayStatus'    Turns on/off the messages generated after checking each
%                   pair of considered in the cross-validation parameters. 
%                   Possible falues are: false and true (default).
%~~~~~~~~~~~~~~~~~~~~~

%%
%-------------------------------------------
%              OUTPUTS:                    -
%-------------------------------------------
% Outputs can be obtained from the structure class object, 'out'.
%~~~~~~~~~~~~~~~~~~~~~
% B                 SpINNEr estimate of B
%~~~~~~~~~~~~~~~~~~~~~
% beta              SpINNEr estimate of beta
%~~~~~~~~~~~~~~~~~~~~~
% bestLambdaN       The best lambda_N found by cross-validation
%~~~~~~~~~~~~~~~~~~~~~
% bestLambdaL       The best lambda_L found by cross-validation
%~~~~~~~~~~~~~~~~~~~~~
% LambsNgrid        Values of lambda_N considered in cross-validation
%~~~~~~~~~~~~~~~~~~~~~
% LambsLgrid        Values of lambda_L considered in cross-validation
%~~~~~~~~~~~~~~~~~~~~~
% logliksCV         The matrix of prediction errors estimated by
%                   cross-validation for each considered pair of 
%                   tuning parameters (lambda_Ns are in rows)
%-------------------------------------------

function out = Logistic_spinnerCV(y, X, AA, W, Params)
%% Tensor of regressor matrices
[p1, p2, p3] = size(AA);
if p3==1
    if mod(p2,p1)~=0
        error('The assumed form is AA = [A1, A2,...,An], with square matrices Ai, hence number of columns in AA should be the multiplicity of the number of rows')
    else
        AA = reshape(AA, [p1, p1, p2/p1]);
    end
end
%% Objects
n = length(y);

%% Cross-validation indices
minElemsN        = floor(n/(Params.kfolds));
remainingElemsN  = n-minElemsN*Params.kfolds;
oneToKfols       = 1:Params.kfolds;
GroupsIdxs       = bsxfun(@times,oneToKfols, ones(minElemsN,1));
GroupsIdxs       = [ GroupsIdxs; [oneToKfols(1:remainingElemsN), zeros(1, Params.kfolds-remainingElemsN)] ]; % if n is not a multiplicity of kfolds
GroupsIdxs       = GroupsIdxs(:);
GroupsIdxs(GroupsIdxs == 0) = [];
s                = rng;
rng('default')  % to get the same results every time
randSmple        = randsample(n,n);
GroupsIdxs       = GroupsIdxs(randSmple);
rng(s);

%% Finding the maximal lambda L
if isempty(Params.gridMaxL)
    clambdaL  = Params.initLambda;
    ValsLambL = zeros(1,1);
    counterr1 = 1;
    stopp     = 0;
    
    % finding lambda_L for which matrix of zeros is obtained
    while stopp == 0
        out = LogisticSpinner(y, X, AA, 0, clambdaL, W, Params);
    %    if norm(out.B, 'fro') < 1e-16
        if sqrt(sum(sum(W.*out.B.^2))) < 1e-16
            stopp = 1;
        end
        ValsLambL(counterr1) = clambdaL;
        clambdaL             = Params.zeroSearchRatio*clambdaL;
        counterr1            = counterr1 + 1;
    end
    
    % initial interval for maximal lambda L
    if isscalar(ValsLambL)
        lamL1 = 0;
        lamL2 = ValsLambL;
    else
        lamL1 = ValsLambL(end-1);
        lamL2 = ValsLambL(end);
    end
    
    % finding narrow interval for maximal lambda L
    stopp      = 0;
    counterr2  = 1;
    ValsLambLmax = zeros(1,1);
    while stopp == 0
        cLamLmaxNew0  = (lamL1 + lamL2)/2;
        outNew0       = LogisticSpinner(y, X, AA, 0, cLamLmaxNew0, W, Params);
        if norm(outNew0.B, 'fro') < 1e-16
            lamL2 = cLamLmaxNew0;
        else
            lamL1 = cLamLmaxNew0;
        end
        ValsLambLmax(counterr2) = lamL2;
        counterr2 = counterr2 + 1;
        if abs(lamL2 - lamL1)/lamL2 < Params.maxLambAcc
            stopp = 1;
        end
    end
else
    lamL2 = Params.gridMaxL;
end

%% Finding the maximal lambda N
if isempty(Params.gridMaxN)
    clambdaN  = Params.initLambda;
    ValsLambN = zeros(1,1);
    counterr1 = 1;
    stopp     = 0;
    
    % finding lambda_N for which matrix of zeros is obtained
    while stopp == 0
        out = LogisticSpinner(y, X, AA, clambdaN, 0, W, Params);
        if norm(out.B, 'fro') < 1e-16
            stopp = 1;
        end
        ValsLambN(counterr1) = clambdaN;
        clambdaN             = Params.zeroSearchRatio*clambdaN;
        counterr1            = counterr1 + 1;
    end
    
    % initial interval for maximal lambda N
    if isscalar(ValsLambN)
        lamN1 = 0;
        lamN2 = ValsLambN;
    else
        lamN1 = ValsLambN(end-1);
        lamN2 = ValsLambN(end);
    end
    
    % finding narrow interval for maximal lambda N
    stopp      = 0;
    counterr2  = 1;
    ValsLambNmax = zeros(1,1);
    while stopp == 0
        cLamLmaxNew0  = (lamN1 + lamN2)/2;
        outNew0       = LogisticSpinner(y, X, AA, cLamLmaxNew0, 0, W, Params);
        if norm(outNew0.B, 'fro') < 1e-16
            lamN2 = cLamLmaxNew0;
        else
            lamN1 = cLamLmaxNew0;
        end
        ValsLambNmax(counterr2) = lamN2;
        counterr2 = counterr2 + 1;
        if abs(lamN2 - lamN1)/lamN2 < Params.maxLambAcc
            stopp = 1;
        end
    end
else
    lamN2 = Params.gridMaxN;
end

%% Final lambdas grids
kL           = Params.expParameterL;
kN           = Params.expParameterN;
seqqL        = (1:(Params.gridLengthL-1))/(Params.gridLengthL-1);
seqqN        = (1:(Params.gridLengthN-1))/(Params.gridLengthN-1);

if Params.addZerosInGrids
    firstValueInGridL = 0;
    firstValueInGridN = 0;
else
    firstValueInGridL = min(1e-6, seqqL(1));
    firstValueInGridN = min(1e-6, seqqN(1));
end

if isequal(Params.searchDensity, 'Old')
    LambsLgrid  = [firstValueInGridL, exp(  ( seqqL * log(lamL2+1)^(1/kL) ).^kL  ) - 1 ];
    LambsNgrid  = [firstValueInGridN, exp(  ( seqqN * log(lamN2+1)^(1/kN) ).^kN  ) - 1 ];
elseif isequal(Params.searchDensity, 'New')
    LambsLgrid  =  [firstValueInGridL, lamL2*(exp(kL.*seqqL)-1 )./(exp(kL) - 1)];
    LambsNgrid  =  [firstValueInGridN, lamN2*(exp(kN.*seqqN)-1 )./(exp(kN) - 1)];
else
    error('Error 451')
end

%% Cross - Validation
logliksCV = zeros(Params.gridLengthN, Params.gridLengthL);
if Params.UseParallel
    if Params.VectorizedParallelization
        logliksCV_3way = zeros(Params.gridLengthN, Params.gridLengthL, Params.kfolds);
        TriplesNumber = Params.gridLengthN * Params.gridLengthL * Params.kfolds;
        parfor idx = 1:TriplesNumber
            [ii, jj, gg] = ind2sub([Params.gridLengthN, Params.gridLengthL, Params.kfolds], idx);
            clambdaN = LambsNgrid(ii);
            clambdaL = LambsLgrid(jj);
            
            % test fold
            testIndices     =  find(GroupsIdxs == gg);
            treningIndices  =  setdiff(1:n , testIndices)';
            AA_trening      =  AA(:,:, treningIndices); %#ok<*PFBNS>
            AA_test         =  AA(:,:, testIndices);
            if ~isempty(X)
                X_trening       =  X(treningIndices,:);
                X_test          =  X(testIndices,:);
            else
                X_trening       =  [];
                X_test          =  0;
            end
            y_trening       =  y(treningIndices);
            y_test          =  y(testIndices);
            out_CV          =  LogisticSpinner(y_trening, X_trening, AA_trening, clambdaN, clambdaL, W, Params);
            AA_test_p       =  permute(AA_test, [3 1 2]);
            eta_test        =  AA_test_p(:,:)*out_CV.B(:) + X_test*out_CV.beta; % correctness confirmed 

            if isequal(Params.scoresFunction, 'LOGLIK')
                logliksCV_3way(idx) =  sum(log(1 + exp(eta_test)) - y_test.*eta_test);
            elseif isequal(Params.scoresFunction, 'NEG_AUC')
                probs = 1./(1 + exp(-eta_test));
                [~,~,~,aucScore] = perfcurve(y_test, probs, 1);
                logliksCV_3way(idx) =  -aucScore;
            elseif isequal(Params.scoresFunction, 'EUCLID')
                probs = 1./(1 + exp(-eta_test));
                logliksCV_3way(idx) =  norm(probs - y_test);
            else
                error('scoresFunction not supported')
            end
        end
        logliksCV = mean(logliksCV_3way,3);
    else
        parfor ii = 1:Params.gridLengthN
        clambdaN = LambsNgrid(ii);
        logliksCVrow = zeros(1, Params.gridLengthL);
            for jj = 1:Params.gridLengthL
                clambdaL  = LambsLgrid(jj);
                normResCV = zeros(Params.kfolds,1);
                for gg = 1:Params.kfolds
                    testIndices     =  find(GroupsIdxs == gg);
                    treningIndices  =  setdiff(1:n , testIndices)';
                    AA_trening      =  AA(:,:, treningIndices); %#ok<*PFBNS>
                    AA_test         =  AA(:,:, testIndices);
                    if ~isempty(X)
                        X_trening       =  X(treningIndices,:);
                        X_test          =  X(testIndices,:);
                    else
                        X_trening       =  [];
                        X_test          =  0;
                    end
                    y_trening       =  y(treningIndices);
                    y_test          =  y(testIndices);
                    out_CV          =  LogisticSpinner(y_trening, X_trening, AA_trening, clambdaN, clambdaL, W, Params);
                    AA_test_p       =  permute(AA_test, [3 1 2]);
                    eta_test        =  AA_test_p(:,:)*out_CV.B(:) + X_test*out_CV.beta; % correctness confirmed 
    
                    if isequal(Params.scoresFunction, 'LOGLIK')
                        normResCV(gg)   =  sum(log(1 + exp(eta_test)) - y_test.*eta_test);
                    elseif isequal(Params.scoresFunction, 'NEG_AUC')
                        probs = 1./(1 + exp(-eta_test));
                        % [xx,yy,~,aa] = perfcurve(y_test, probs, 1);
                        [~,~,~,aucScore] = perfcurve(y_test, probs, 1);
                        normResCV(gg)   =  -aucScore;
                    elseif isequal(Params.scoresFunction, 'EUCLID')
                        probs = 1./(1 + exp(-eta_test));
                        normResCV(gg) =  norm(probs - y_test);
                    else
                        error('scoresFunction not supported')
                    end
                end
                logliksCVrow(jj) =  sum(normResCV)/Params.kfolds;
                if Params.displayStatus
                    disp(strcat(['finished:  ', num2str(ii), endText(ii),' ', 'value from lambdaN grid,  ', num2str(jj), endText(jj), ' value from lambdaL grid']))
                end
            end
            logliksCV(ii, :) = logliksCVrow;
        end
    end
else
    for ii = 1:Params.gridLengthN
        clambdaN = LambsNgrid(ii);
        for jj = 1:Params.gridLengthL
            clambdaL  = LambsLgrid(jj);
            normResCV = zeros(Params.kfolds,1);
            for gg = 1:Params.kfolds
                testIndices     =  find(GroupsIdxs == gg);
                treningIndices  =  setdiff(1:n , testIndices)';
                AA_trening      =  AA(:,:, treningIndices); %#ok<*PFBNS>
                AA_test         =  AA(:,:, testIndices);
                if ~isempty(X)
                    X_trening       =  X(treningIndices,:);
                    X_test          =  X(testIndices,:);
                else
                    X_trening       =  [];
                    X_test          =  0;
                end
                y_trening       =  y(treningIndices);
                y_test          =  y(testIndices);
                out_CV          =  LogisticSpinner(y_trening, X_trening, AA_trening, clambdaN, clambdaL, W, Params);
                AA_test_p       =  permute(AA_test, [3 1 2]);
                eta_test        =  AA_test_p(:,:)*out_CV.B(:) + X_test*out_CV.beta;
                
                if isequal(Params.scoresFunction, 'LOGLIK')
                    normResCV(gg)   =  sum(log(1 + exp(eta_test)) - y_test.*eta_test);
                elseif isequal(Params.scoresFunction, 'NEG_AUC')
                    probs = 1./(1 + exp(-eta_test));
                    % [xx,yy,~,aa] = perfcurve(y_test, probs, 1);
                    [~,~,~,aucScore] = perfcurve(y_test, probs, 1);
                    normResCV(gg)   =  -aucScore;
                elseif isequal(Params.scoresFunction, 'EUCLID')
                    probs = 1./(1 + exp(-eta_test));
                    normResCV(gg) =  norm(probs - y_test);
                else
                    error('scoresFunction not supported')
                end
            end        
            logliksCV(ii, jj) =  sum(normResCV)/Params.kfolds;
            if Params.displayStatus
                disp(strcat(['finished:  ', num2str(ii), endText(ii),' ', 'value from lambdaN grid,  ', num2str(jj), endText(jj), ' value from lambdaL grid']))
            end
        end
    end
end

%% Optimal lambdas
[MM, II]    = min(logliksCV);
Xindex      = find(MM == min(MM), 1);
Yindex      = II(Xindex);
bestLambdaN = LambsNgrid(Yindex);
bestLambdaL = LambsLgrid(Xindex);

%% Final estimate
outFinal    =  LogisticSpinner(y, X, AA, bestLambdaN, bestLambdaL, W, Params);

%% Output
out              = struct;
out.LambsLgrid   = LambsLgrid;
out.LambsNgrid   = LambsNgrid;
out.logliksCV    = logliksCV;
out.B            = outFinal.B;
out.beta         = outFinal.beta;
out.bestLambdaN  = bestLambdaN;
out.bestLambdaL  = bestLambdaL;
out.GroupsIdxs   = GroupsIdxs;
out.optVal       = outFinal.optVal;
out.optVal_lamNucl  = outFinal.optVal_lamNucl;
out.optVal_loglik   = outFinal.optVal_loglik;
out.optVal_lamLasso = outFinal.optVal_lamLasso;

% EstimatedProbs and AUC
AA2       =  permute(AA, [3 1 2]);
eta_final  =  AA2(:,:)*out.B(:) + X*out.beta;
probs_final = 1./(1 + exp(-eta_final));
% [xx,yy,~,aa] = perfcurve(y, probs_final, 1);
[~,~,~,AUC] = perfcurve(y, probs_final, 1);
out.AUC = AUC;
out.Probs = probs_final;
out.Euclid = norm(y - probs_final);

end