%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB programme file for the toolkit for           %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from the Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: GMA, 01/2024                           %%%
% Last updated by GMA 06/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function is not part of the original replication directory      %%%
%%% This function merges a block-level outcome to a shapefile of         %%%
%%% statistical blocks, generates a map, and saves it at a desired       %%%
%%% destination                                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program uses the following inputs:
% shapefile - the block shapefile to which the outputs are being merged
% variable  - the outcome you wish to map  
% name      - the title you wish to appear on the map    
% folder    - the location where you wish to save the map
% filename  - the filename of the map to be generated

% The program produces the following outputs:
% RESULT - a placeholder used to display confirmation of successful completion
% A map with the desired filename will be saved in the selected folder in PNG format

function RESULT = MAPIT(shapefile, variable, name, folder, filename)
    % Clear any pre-existing figure
    clf;
    close all;

    % Define the number of colors you want in your colormap
    numColors = 256; % Typical colormap length

    % Create a matrix of RGB values that transition from yellow to red
    customColormap = [ones(numColors, 1) linspace(1, 0, numColors)' zeros(numColors, 1)];

    % Read shapefiles
    GREEN = shaperead('../shapefile/BerlinGreen');  % Green spaces shapefile
    WATER = shaperead('../shapefile/BerlinWater');  % Water spaces shapefile
    SHAPE = shaperead(shapefile);                   % User-specified shapefile
    BOUNDARIES = shaperead('../shapefile/Bezirke23'); % Bezirke23 boundaries shapefile

    % Calculate Jenks breaks
    numClasses = 10;  % Number of classes
    breaks = jenks(variable, numClasses);

    % Assign the outcome variable and categorize based on Jenks breaks
    for i = 1:numel(SHAPE)
        if isnan(variable(i)) || variable(i) == 0
            SHAPE(i).Category = 0; % Category for 0 or NaN values
        else
            SHAPE(i).SHOW = variable(i);
            for d = 1:length(breaks)-1
                if variable(i) <= breaks(d+1)
                    SHAPE(i).Category = d;
                    break;
                end
            end
        end
    end

    % Define the symbol specification for the map with categories
    colorIndices = round(linspace(1, numColors, numClasses)); % Choose colors from the colormap
    symbolSpec = makesymbolspec('Polygon', ...
        {'Category', 0, 'FaceColor', [1 1 1]}, ... % White for 0 or NaN values
        {'Category', 1, 'FaceColor', customColormap(colorIndices(1), :)}, ...
        {'Category', 2, 'FaceColor', customColormap(colorIndices(2), :)}, ...
        {'Category', 3, 'FaceColor', customColormap(colorIndices(3), :)}, ...
        {'Category', 4, 'FaceColor', customColormap(colorIndices(4), :)}, ...
        {'Category', 5, 'FaceColor', customColormap(colorIndices(5), :)}, ...
        {'Category', 6, 'FaceColor', customColormap(colorIndices(6), :)}, ...
        {'Category', 7, 'FaceColor', customColormap(colorIndices(7), :)}, ...
        {'Category', 8, 'FaceColor', customColormap(colorIndices(8), :)}, ...
        {'Category', 9, 'FaceColor', customColormap(colorIndices(9), :)}, ...
        {'Category', 10, 'FaceColor', customColormap(colorIndices(10), :)}, ...
        {'Default', 'EdgeColor', 'none'});

    % Generate the map figure using the specified shapefile
    figure; % Open a new figure window
    mapshow(SHAPE, 'SymbolSpec', symbolSpec);
    hold on;
    mapshow(WATER, 'FaceColor', [0.678 0.847 0.902], 'EdgeColor', 'none'); % Add water shape
    mapshow(GREEN, 'FaceColor', [0.133 0.545 0.133], 'EdgeColor', 'none'); % Add green shape
    mapshow(BOUNDARIES, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 0.25); % Add Bezirke23 boundaries with even thinner lines

    % Create a colorbar that matches the colors in the map
    colormap([[1 1 1]; customColormap(colorIndices, :)]); % Include white for the first category
    
    % Create tick labels for colorbar showing ranges
    tickLabels = {'0 or NaN'};
    for k = 1:numClasses
        tickLabels{end+1} = sprintf('%.2f - %.2f', breaks(k), breaks(k+1));
    end

    % Create the colorbar
    hcb = colorbar('Ticks', (0.5:1:(numClasses+0.5)), 'TickLabels', tickLabels);
    caxis([0 numClasses+1]); % Set the limits of the colorbar to match the number of classes plus the extra category
    set(get(hcb, 'Title'), 'String', 'Jenks Breaks');

    % Add title to the map
    title(name);

    % Generate the file path
    filepath = fullfile(folder, [filename, '.png']);
    
    % Save the figure with higher resolution
    print(gcf, filepath, '-dpng', '-r300'); % Save as PNG with 300 DPI

    % Confirm successful execution
    RESULT = 'DONE';
end

% Jenks Natural Breaks Algorithm Implementation
function breaks = jenks(data, numClasses)
    data = data(~isnan(data) & data ~= 0); % Exclude NaN and 0 values
    data = sort(data(:));
    numData = length(data);
    
    % Initialize matrices
    lowerClassLimits = ones(numData, numClasses);
    varianceCombinations = inf(numData, numClasses);
    
    % Calculate variance for each potential class
    for i = 2:numData
        sum1 = 0;
        sum2 = 0;
        w = 0;
        
        for j = 1:i
            k = i - j + 1;
            sum1 = sum1 + data(k);
            sum2 = sum2 + data(k) * data(k);
            w = w + 1;
            variance = sum2 - (sum1 * sum1) / w;
            
            if k ~= 1
                for m = 2:numClasses
                    if varianceCombinations(i, m) >= (variance + varianceCombinations(k - 1, m - 1))
                        lowerClassLimits(i, m) = k;
                        varianceCombinations(i, m) = variance + varianceCombinations(k - 1, m - 1);
                    end
                end
            end
        end
        lowerClassLimits(i, 1) = 1;
        varianceCombinations(i, 1) = variance;
    end
    
    % Find the k-class
    breaks = zeros(1, numClasses);
    k = numData;
    for j = numClasses:-1:1
        breaks(j) = data(k);
        k = lowerClassLimits(k, j) - 1;
    end
    breaks = [data(1), breaks];
end
