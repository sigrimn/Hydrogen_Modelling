function modelCoarse = makeCoarseModel2D(modelRef, dimensions)

map = modelRef.G.cells.indexMap;

% Nummerer indexMap unikt
modelRef.G.cells.indexMap = (1:numel(modelRef.G.cells.indexMap))';

blockIx = partitionUI(modelRef.G, dimensions);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);
% Perform a simple upscaling to obtain a coarse model
modelRef.G.cells.indexMap = map;
modelCoarse = upscaleModelTPFA(modelRef, blockIx);
modelCoarse.AutoDiffBackend = AutoDiffBackend();
% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).

% Antall rader du ønsker
numRows = 3;

pts = modelCoarse.fluid.krPts;

scaling = {'SGL',   pts.g(1), 'SGCR', pts.g(2), 'SGU', pts.g(3), ...
            'SOGCR', pts.og(2), 'KRG',  pts.g(4), 'KRO', pts.og(4)};
modelCoarse = imposeRelpermScaling(modelCoarse, scaling{:});
modelCoarse.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

% % Utvid g og og til flere rader
modelCoarse.fluid.krPts.g = repmat(pts.g, numRows, 1);
modelCoarse.fluid.krPts.og = repmat(pts.og, numRows, 1);

% Manually setting upper left to obtain symmetry
modelCoarse.rock.regions.saturation((6 * 11 + 2)) = 3;

% Validating 
modelCoarse = modelCoarse.validateModel();