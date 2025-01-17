function modelCoarse = makeCoarseModel3D(modelRef, dimensions)

blockIx = partitionUI(modelRef.G, dimensions);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);
% Perform a simple upscaling to obtain a coarse model
modelCoarse = upscaleModelTPFA(modelRef, blockIx);
modelCoarse.AutoDiffBackend = AutoDiffBackend();
% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).

%pts = modelCoarse.fluid.krPts;
%scaling = {'SGL',   pts.g(1), 'SGCR', pts.g(2), 'SGU', pts.g(3), ...
%            'SOGCR', pts.og(2), 'KRG',  pts.g(4), 'KRO', pts.og(4)};
%modelCoarse = imposeRelpermScaling(modelCoarse, scaling{:});
modelCoarse.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

% Validating 
modelCoarse = modelCoarse.validateModel();