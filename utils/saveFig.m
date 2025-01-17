function varargout = saveFig(figure, filename, outputDir,varargin)

figure.Color = [1, 1, 1];
% Define the file path to a different folder
exportgraphics(figure, fullfile(outputDir, filename), 'ContentType', 'vector');

