classdef spinner < handle
    properties
        y
        X
        A % matrix of a form [A1,A2,...]
        AA % 3- way tensor
        W

        % sizes
        n
        p
        d
        modelFitted = false
        
        % available
        Estimate
        Hclust = struct('POS', struct('Similarity', [], 'Distance', [], 'Bsign', [], 'Hclusters', [], 'LeafOrder', [], 'InvLeafOrder', [], 'Bthresh1', [], 'Breord', [], 'Bmasked', [], 'Bgray', [], 'Bbin', [], 'Bclusts', []), ...
            'NEG', struct('Similarity', [], 'Distance', [], 'Bsign', [], 'Hclusters', [],  'LeafOrder', [], 'InvLeafOrder', [], 'Bthresh1', [], 'Breord', [], 'Bmasked', [], 'Bgray', [], 'Bbin', [], 'Bclusts', []),...
            'Clusters', struct());
        Clusters = struct('POS', struct(), 'NEG', struct());
        History = {}
    end

    methods
        %------------------------------------------------------------------
        % Initialization
        %------------------------------------------------------------------
        function obj = spinner(y, X, A, W)
            if nargin == 0
                % show all the list of all possible arguments ( and the corresponding default values)
                defaultSettings = ParseArgumentsFit({});
                disp(defaultSettings)
                return;
            elseif nargin < 4
                StandardWeightsMatrix = true; % 'Zeros on diagonal, ones on off-diagonal'
            else
                StandardWeightsMatrix = false;
            end

            % Assign the required arguments
            obj.y = y;
            obj.X = X;
            obj.d = size(X,2);

            % Tensor of regressor matrices
            obj.convertRegressorMatrices(A);

            if StandardWeightsMatrix
                obj.W = ones(obj.p, obj.p) - eye(obj.p);
            else
                obj.W = W;
            end
        end

        %------------------------------------------------------------------
        function convertRegressorMatrices(obj, A)
            [p1, p2, p3] = size(A);
            if p3==1
                if mod(p2,p1)~=0
                    error('The assumed form is AA = [A1, A2,...,An], with square matrices Ai, hence number of columns in AA should be the multiplicity of the number of rows')
                else
                    obj.A = A;
                    obj.AA = reshape(A, [p1, p1, p2/p1]);
                end
            else
                if p1~=p2
                    error('The matrices Ai are assumed to be square')
                else
                    obj.A = reshape(A, [p1*p2, p3])';
                    obj.AA = A;
                end
            end

            % sizes
            obj.p = p1;
            obj.n = p3;
        end

        %------------------------------------------------------------------
        function out = fit(obj, varargin)
            % Collect input arguments
            [ParamsFit, ~] = ParseArgumentsFit(varargin);

            %--------------------------------------------------------------
            %                  PROVIDED TUNING PARAMETERS
            %--------------------------------------------------------------

            if ~isempty(ParamsFit.LambdaL) || ~isempty(ParamsFit.LambdaN)
                if isempty(ParamsFit.LambdaL) || isempty(ParamsFit.LambdaN)
                    if ~isempty(ParamsFit.LambdaL)
                        % CV for single nuclear norm penalty
                        if isequal(ParamsFit.Family, 'Gaussian')
                            if isequal(ParamsFit.Method, 'CV')
                                % SPINNER: Find parameter by Cross-Validation
                                error('CV Not implemented')
                                % out = spinnerCVnuclear(obj.y, obj.X, obj.AA, obj.W, ParamsFit);  <---- just a placeholder
                            else
                                % SPINNER: Find parameter by Bayesian
                                error('Not implemented')
                            end
                        elseif isequal(ParamsFit.Family, 'Binomial')
                            assert(all(or(obj.y==0, obj.y == 1)), 'Response vector, y, should have only zeros and ones coefficients if Family == "Binomial"')
                            if isequal(ParamsFit.Method, 'CV')
                                % LOGISTIC SPINNER: Find parameters by Cross-Validation
                                out = Logistic_spinnerCVnuclear(obj.y, obj.X, obj.AA, obj.W, ParamsFit);
                            else
                                error('Not implemented')
                            end
                        else
                            error('Not implemented')
                        end
                    else
                        % CV for single lasso norm penalty
                        if isequal(ParamsFit.Family, 'Gaussian')
                            if isequal(ParamsFit.Method, 'CV')
                                % SPINNER: Find parameter by Cross-Validation
                                % out = spinnerCVlasso(obj.y, obj.X, obj.AA, obj.W, ParamsFit);
                            else
                                % SPINNER: Find parameter by Bayesian
                                error('Not implemented')
                            end
                        elseif isequal(ParamsFit.Family, 'Binomial')
                            assert(all(or(obj.y==0, obj.y == 1)), 'Response vector, y, should have only zeros and ones coefficients if Family == "Binomial"')
                            if isequal(ParamsFit.Method, 'CV')
                                % LOGISTIC SPINNER: Find parameter by Cross-Validation
                                out = Logistic_spinnerCVlasso(obj.y, obj.X, obj.AA, obj.W, ParamsFit);
                            else
                                error('Not implemented')
                            end
                        else
                            error('Not implemented')
                        end
                    end
                    % error('[not yet implemented] both or none parameters must be provided')
                else
                    % two parameters provided
                    if isequal(ParamsFit.Family, 'Gaussian')
                        out = spinnerRun(obj.y, obj.X, obj.AA, ParamsFit.LambdaN, ParamsFit.LambdaL, obj.W, ParamsFit);
                    elseif isequal(ParamsFit.Family, 'Binomial')
                        assert(all(or(obj.y==0, obj.y == 1)), 'Response vector, y, should have only zeros and ones coefficients if Family == "Binomial"')
                        out = LogisticSpinner(obj.y, obj.X, obj.AA, ParamsFit.LambdaN, ParamsFit.LambdaL, obj.W, ParamsFit);
                    else
                        error('Not implemented')
                    end
                end
            else

            %--------------------------------------------------------------
            %                 NO TUNING PARAMETERS
            %--------------------------------------------------------------

                if isequal(ParamsFit.Family, 'Gaussian')
                    if isequal(ParamsFit.Method, 'CV')
                        % SPINNER: Find parameters by Cross-Validation
                        out = spinnerCV(obj.y, obj.X, obj.AA, obj.W, ParamsFit);
                    else
                        % SPINNER: Find parameters by Bayesian
                        error('Not implemented')
                    end
                elseif isequal(ParamsFit.Family, 'Binomial')
                    assert(all(or(obj.y==0, obj.y == 1)), 'Response vector, y, should have only zeros and ones coefficients if Family == "Binomial"')
                    if isequal(ParamsFit.Method, 'CV')
                        % LOGISTIC SPINNER: Find parameters by Cross-Validation
                        out = Logistic_spinnerCV(obj.y, obj.X, obj.AA, obj.W, ParamsFit);
                    elseif isequal(ParamsFit.Method, 'CV_COORD')
                        % %----------------------------------------------------------------------
                        % % CV for nuclear norm
                        % ParamsFit_N = ParamsFit; ParamsFit_N.LambdaL = 0;
                        % out_N = Logistic_spinnerCVnuclear(obj.y, obj.X, obj.AA, obj.W, ParamsFit_N);
                        % % CV for lasso norm
                        % ParamsFit_L = ParamsFit; ParamsFit_L.LambdaN = out_N.bestLambdaN;
                        % out_L = Logistic_spinnerCVlasso(obj.y, obj.X, obj.AA, obj.W, ParamsFit_L);
                        % % Final estimate
                        % ParamsFit.LambdaN = out_N.bestLambdaN; ParamsFit.LambdaL = out_L.bestLambdaL;
                        % assert(all(or(obj.y==0, obj.y == 1)), ['Response vector, y, ' ...
                        %     'should have only zeros and ones coefficients if Family == "Binomial"'])
                        % out = LogisticSpinner(obj.y, obj.X, obj.AA, ParamsFit.LambdaN, ParamsFit.LambdaL, obj.W, ParamsFit);
                        % 
                        %----------------------------------------------------------------------
                        % CV for lasso norm
                        ParamsFit_L = ParamsFit; ParamsFit_L.LambdaN = 0;
                        out_L = Logistic_spinnerCVlasso(obj.y, obj.X, obj.AA, obj.W, ParamsFit_L);
                        % CV for nuclear norm
                        ParamsFit_N = ParamsFit; ParamsFit_N.LambdaL = out_L.bestLambdaL;
                        out_N = Logistic_spinnerCVnuclear(obj.y, obj.X, obj.AA, obj.W, ParamsFit_N);
                        % Final estimate
                        ParamsFit.LambdaN = out_N.bestLambdaN; ParamsFit.LambdaL = out_L.bestLambdaL;
                        assert(all(or(obj.y==0, obj.y == 1)), ['Response vector, y, ' ...
                            'should have only zeros and ones coefficients if Family == "Binomial"'])
                        out = LogisticSpinner(obj.y, obj.X, obj.AA, ParamsFit.LambdaN, ParamsFit.LambdaL, obj.W, ParamsFit);
                        %----------------------------------------------------------------------
                    else
                        error('Not implemented')
                    end
                else
                    error('Not implemented')
                end
            end
            obj.Estimate = out;
            obj.History = cat(1, obj.History, out);
            obj.modelFitted = true;
        end

        %------------------------------------------------------------------
        function identifyClusters2(obj, clustNum)
            if ~obj.modelFitted
                error('To plot the estimate, the model must be fitted first by using "fit" method.')
            end
            % k = 5
            [V, ~] = eigs(obj.Estimate.B, clustNum+1);  
            clusterLabels = kmeans(V, clustNum+1);
            
            %Louvain community detection algorithm
            % https://github.com/GenLouvain/GenLouvain
            % https://sites.google.com/site/bctnet/network-visualization

        end

        %------------------------------------------------------------------
        function identifyClusters(obj, expected_max_num_blocks, expected_max_blocks_size)
            % Recover block-diagonal structure from a symmetric matrix A
            % threshold: tolerance for small off-block entries (default = 0.1)

            if ~obj.modelFitted
                error('To plot the estimate, the model must be fitted first by using "fit" method.')
            end
            if nargin < 2
                expected_max_num_blocks = 10;
            end
            if nargin < 3
                expected_max_blocks_size = Inf;
            end

            %--------------------------------------------------------------
            %              PARAMETERS
            %--------------------------------------------------------------
            THRESH_PROP = 0.3;
            THRESH_PROP2 = 0.3;
            NODES_PERM = 'OPT_LEAF_ORDER';
            % NODES_PERM = 'FROM_CLUSTERS';
            LINKAGE_TYPE = 'average';
            % LINKAGE_TYPE = 'complete';
            % LINKAGE_TYPE = 'ward';
            FILL_HOLES = true;


            assert(ismember(NODES_PERM, {'OPT_LEAF_ORDER', 'FROM_CLUSTERS'}))
            assert(ismember(LINKAGE_TYPE, {'average', 'complete', 'ward'}))

            %--------------------------------------------------------------
            % create bandwidth filter
            [row, col] = ndgrid(1:obj.p, 1:obj.p);
            bandWidthMask = abs(row - col) < expected_max_blocks_size;

            for effect = {'POS', 'NEG'}
                ef = effect{1};

                % Compute similarity matrix
                cB = obj.Estimate.B;
                if isequal(effect{1}, 'POS')
                    cB(obj.Estimate.B<0) = 0;
                elseif isequal(effect{1}, 'NEG')
                    cB(obj.Estimate.B>0) = 0;
                else
                    error('Effect not recognized')
                end
                obj.Hclust.(ef).Bsign = cB;

                % Thresholding
                cB = abs(cB);
                maxEntry = max(max(cB));
                FiltThresh = THRESH_PROP*maxEntry;
                cB(cB < FiltThresh) = 0;
                obj.Hclust.(ef).Bthresh1 = cB;
    
                % Similarity matrix
                similarity = cB;
                similarity = similarity - min(min(similarity));
                similarity(similarity>1) = 1;
                similarity = similarity/max(max(similarity));
                obj.Hclust.(ef).Similarity = similarity;

                % Distance matrix
                distance = 1 - similarity;
                distance = distance - diag(diag(distance));
                obj.Hclust.(ef).Distance = distance;
            
                % Hierarchical clustering: leaf order
                Y = squareform(distance);
                Z = linkage(Y, LINKAGE_TYPE);
                if isequal(NODES_PERM, 'OPT_LEAF_ORDER')
                    leafOrder = optimalleaforder(Z, Y);
                elseif isequal(NODES_PERM, 'FROM_CLUSTERS')
                    clusters = cluster(Z, 'maxclust', expected_max_num_blocks);
                    [~,leafOrder] = sort(clusters);
                else
                    error('NODES_PERM should be equal to "OPT_LEAF_ORDER" or "FROM_CLUSTERS"')
                end
                
                obj.Hclust.(ef).LeafOrder = leafOrder;
                invLeafOrder = zeros(size(leafOrder));
                invLeafOrder(leafOrder) = 1:length(leafOrder);
                obj.Hclust.(ef).InvLeafOrder = invLeafOrder;
                if isequal(NODES_PERM, 'FROM_CLUSTERS')
                    obj.Hclust.(ef).Hclusters = clusters;
                end

                % Reordered estimate
                cBreord = cB(leafOrder, leafOrder);
                obj.Hclust.(ef).Breord = cBreord;

                % Masking with widthBand filter
                cBreordMasked = cBreord;
                cBreordMasked(~bandWidthMask) = 0;
                obj.Hclust.(ef).Bmasked = cBreordMasked;

                % Converting to grayshades entries
                cBreordMaskedGray = cBreordMasked-(diag(diag(cBreordMasked)));
                cBreordMaskedGray = cBreordMaskedGray/max(max(cBreordMaskedGray));

                % Second filtering
                cBreordMaskedGray(cBreordMaskedGray<THRESH_PROP2) = 0;

                % Fill holes (optional)
                if FILL_HOLES
                    cBreordMaskedGray = imfill(cBreordMaskedGray);
                end
                obj.Hclust.(ef).Bgray = cBreordMaskedGray;

                % Converting to binary entries
                cBreordMaskedBin = cBreordMaskedGray > 1e-8;
                obj.Hclust.(ef).Bbin = cBreordMaskedBin;
                
                if false
                % Finding structures
                connComp = bwconncomp(cBreordMaskedBin, 4);
                RowsIdxsOrdered = [];
                ColsIdxsOrdered = [];
                RowsIdxs = [];
                ColsIdxs = [];
                IsDiagonal = [];
                LinearIdxsOrdered = [];
                LinearIdxs = [];
                ID = [];
                ccounter = 1;
                for ss = 1:length(connComp.PixelIdxList)
                    if length(connComp.PixelIdxList{ss}) > 1
                        [crows, ccols] = ind2sub([obj.p,obj.p], connComp.PixelIdxList{ss});
                        crowsOrig = leafOrder(crows);
                        ccolsOrig = leafOrder(ccols);
                        RowsIdxsOrdered = [RowsIdxsOrdered;crows]; %#ok<AGROW>
                        ColsIdxsOrdered = [ColsIdxsOrdered;ccols]; %#ok<AGROW>
                        RowsIdxs = [RowsIdxs;crowsOrig]; %#ok<AGROW>
                        ColsIdxs = [ColsIdxs;ccolsOrig]; %#ok<AGROW>
                        if any(crows==ccols)
                            IsDiagonal = [IsDiagonal; ones(length(crows),1)]; %#ok<AGROW>
                        else
                            IsDiagonal = [IsDiagonal; zeros(length(crows),1)]; %#ok<AGROW>
                        end
                        linearIndex = sub2ind([obj.p,obj.p], crowsOrig, ccolsOrig);
                        LinearIdxsOrdered = [LinearIdxsOrdered; connComp.PixelIdxList{ss}]; %#ok<AGROW>
                        LinearIdxs = [LinearIdxs; linearIndex]; %#ok<AGROW>
                        ID = [ID;ccounter*ones(length(crows),1)]; %#ok<AGROW>
                        ccounter = ccounter + 1;
                    end
                end
                structTable = table(ID, RowsIdxsOrdered, ColsIdxsOrdered, ...
                    LinearIdxsOrdered, RowsIdxs, ColsIdxs, LinearIdxs, IsDiagonal);
                obj.Clusters.(effect{1}).StructuresSummary = structTable;

                % identify symmetric off-diagonal counterparts
                ID2 = zeros(length(ID),1);
                remIdxsSet = 1:max(ID);
                stopp = 0;
                cNewId = 1;
                while stopp == 0
                    cidx = remIdxsSet(1);
                    cidx2 = obj.identifySymmetricID(effect{1}, cidx);
                    ID2(or(ID==cidx, ID==cidx2)) = cNewId;
                    remIdxsSet = setdiff(remIdxsSet, [cidx, cidx2]);
                    cNewId = cNewId + 1;
                    if isempty(remIdxsSet)
                        stopp = 1;
                    end
                end
                obj.Clusters.(effect{1}).StructuresSummary.ID2 = ID2;
                structMatrix = sparse(RowsIdxs, ColsIdxs, ID2, obj.p, obj.p);
                obj.Clusters.(effect{1}).StructuresMatrix = structMatrix;
                end
            end
            % obj.recoverBlocks()
        end

        %------------------------------------------------------------------
        function [recoveredMatrix, permutation] = identifyClusters3(obj, expected_num_blocks)
            % Recover block-diagonal structure from a symmetric matrix A
            % threshold: tolerance for small off-block entries (default = 0.1)
        
            if nargin < 2
                expected_num_blocks = 10;
            end
            % ----------------------------------------------------------
            %                   POSITIVE ENTRIES
            % ----------------------------------------------------------
            
            % Compute similarity matrix
            cB = obj.Estimate.B; 
            cB(obj.Estimate.B<0) = 0;
            similarity = abs(cB);

            % Scale similarity matrix and generate distance matrix
            similarity = similarity - min(min(similarity));
            similarity(similarity>1) = 1;
            similarity = similarity/max(max(similarity));
            distance = 1 - similarity;
            distance = distance - diag(diag(distance));
        
            % DBSCAN
            blockSize = 10;
            MINPTS = 5;
            EPSILON = blockSize;
            Y = squareform(distance); % Convert to condensed form
            % IDX = dbscan(Y, EPSILON, MINPTS, 'DISTANCE', 'PRECOMPUTED') ;

            %  DBSCAN2
            % Threshold to ignore noise
            threshold = 0.01;
            [row, col] = find(abs(cB) > threshold);
            points = [row, col];
            epsilon = 0.01; % max distance between points in a cluster
            minPts = 2;          % minimum number of points to form a cluster
            labels = dbscan(points, epsilon, minPts);
            unique(labels)

        end

        %------------------------------------------------------------------
        function recoverBlocks(obj)
            POStable = obj.Clusters.POS.StructuresSummary;
            NEGtable = obj.Clusters.NEG.StructuresSummary;
            
            % restrict to the selected columns
            selColNames = {'ID2', 'RowsIdxsOrdered','ColsIdxsOrdered','RowsIdxs','ColsIdxs', 'IsDiagonal'};
            POSsubtab = POStable(:, selColNames);
            POSsubtab.Positive = ones(height(POSsubtab),1);
            NEGsubtab = NEGtable(:, selColNames);
            NEGsubtab2 = NEGsubtab;
            NEGsubtab2.ID2 = NEGsubtab2.ID2 + max(POSsubtab.ID2);
            NEGsubtab2.Positive = zeros(height(NEGsubtab2),1);
            
            % generate total summary tables for identified clusters
            totTable = [POSsubtab; NEGsubtab2];
            
            Nall = max(totTable.ID2);
            NumPixels = zeros(Nall,1);
            IsDiag = zeros(Nall,1);
            BlocksIdxs = cell(Nall,1);
            BlocksSizes = zeros(Nall,1);
            for ii = 1:Nall
                cidxs = totTable.ID2 == ii;
                NumPixels(ii) = sum(cidxs);
                IsDiag(ii) = unique(totTable.IsDiagonal(cidxs));
                BlocksIdxs{ii} = unique([totTable.RowsIdxs(cidxs);totTable.ColsIdxs(cidxs)], 'stable');
                BlocksSizes(ii) = length(BlocksIdxs{ii});
            end

            % create Indicator Matrix (IM) for blocks
            IMrows = vertcat(BlocksIdxs{:});
            IMcols = repelem(1:Nall, BlocksSizes)';
            IMvals = ones(length(IMcols),1);
            IM = sparse(IMrows, IMcols, IMvals);

            % ------------------------------
            % TEMPORARY
            % -------------------------------
            % cidx = totTable.ID2 == 15;
            % 
            % aa=[BlocksIdxs{2};BlocksIdxs{16}];
            % 
            % aa=[BlocksIdxs{2};BlocksIdxs{5};BlocksIdxs{8};BlocksIdxs{16}];
            % 
            % spinnerHeatmap(obj.Estimate.B(aa,aa))
            % 
            % % 
            % % [totTable.RowsIdxs(cidx);totTable.ColsIdxs(cidx)]
            % 
            % uniqueValues = unique([totTable.RowsIdxs(cidx);totTable.ColsIdxs(cidx)], 'stable');
            % spinnerHeatmap(obj.Estimate.B(uniqueValues,uniqueValues))
        end

        %------------------------------------------------------------------
        function ID = identifySymmetricID(obj, pos_neg, id)
            if ~ismember(pos_neg, {'POS', 'NEG'})
                error('"pos_neg" must be either "POS" or "NEG"')
            end
            
            structTable = obj.Clusters.(pos_neg).StructuresSummary;
            IDnum = max(structTable.ID);

            if id > IDnum
                error('The selected "id" can not be larger than the number of all identified structures')
            end

            structureIdxs = structTable.ID == id;
            if ~unique(structTable.IsDiagonal(structureIdxs)) == 0
                % disp('The selected structure is on-diagonal - no need to find the symmetric counterpart')
                ID = id;
                return
            end

            referenceRows = sort(structTable.RowsIdxsOrdered(structureIdxs));
            referenceCols = sort(structTable.ColsIdxsOrdered(structureIdxs));

            for ii = 1:IDnum
                cIdxs = structTable.ID == ii;
                cRows = sort(structTable.RowsIdxsOrdered(cIdxs));
                cCols = sort(structTable.ColsIdxsOrdered(cIdxs));
                if and(isequal(cRows, referenceCols), isequal(cCols, referenceRows))
                    ID = ii;
                    return
                end
            end
            error('Symmetric counterpart was not found')

        end

        %------------------------------------------------------------------
        function plotEstimate(obj, histCounter)
            if nargin == 1
                histCounter = length(obj.History);
            end

            if ~obj.modelFitted
                error('To plot the estimate, the model must be fitted first by using "fit" method.')
            end

            % get estimate and create its heatmap
            b = obj.History{histCounter}.B;
            b = double(b);      
            bottom  = min(min(b));  % fixed option for now
            top = max(max(b));  % fixed option for now
            tit =''; % fixed option for now
            margins = [0.05, 0.05];
            figure('name','Estimates');
            axes('Position', [margins(2), margins(1), (1-2*margins(2)), (1-2*margins(1))]);
            imagesc(b);
            colormap(spinnerColormap(bottom ,top));
            title(tit);
            axis equal;
            axis tight;
            clim manual
            clim([bottom  top]);
            colorbar
        end

        %------------------------------------------------------------------
        function updateWeights(obj, W)
            obj.W = W;
        end

        %------------------------------------------------------------------
        function subStruct = getSubstruct(~, S, fields)
            Fields = fieldnames(S);
            foundFieldsIdxs = ismember(fields,Fields);
            if ~all(foundFieldsIdxs)
                missingFields = join(fields(~foundFieldsIdxs),', ');
                error(['The following fields names were not identified in struct: ', missingFields{1}])
            end
            
            nFields = length(fields);
            for ii = 1:nFields
                if ismember(fields{ii}, Fields)
                    subStruct.(fields{ii}) = S.(fields{ii});
                end
            end
        end
        
    end
end
