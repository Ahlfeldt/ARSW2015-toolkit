%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB programme file for the toolkit for           %%%
%%% Ahlfeldt, Redding, Sturm, Wolf (2015)               %%%
%%% Economics of density: Evidence from teh Berlin Wall %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version: GMA, 01/2024                           %%%
% Last updated by GMA 03/2024                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function is not part of the orginal replication directory      %%%
%%% This function merges a block-level outcome to a shapefile of        %%%
%%% statistical blocks, generates a map, and saves it at a desired      %%%
%%% destination                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This below program uses the following inputs
    % shapefile is the the block shapefile to which the outputs are beign merged
            % This has to be the Berlin4matlab shapefile. It has been
            % resorted so to use the same indexing as the arrays in this
            % workspace
        % variable is the outcome you wish to map  
        % name is the title you wish appear on the map    
        % folder is the location where you wish to save the map
        % filename is the filename of the map to be generated
% The below program produces the following outputs
        % RESULT is a placeholder used to display confirmation of
            % successful cmpletion
        % A map with the desired filename will be saved in the selected
            % folder in png format
function RESULT = MAPIT(shapefile,variable,name,folder,filename)
clf;
close all;
% Define the number of colors you want in your colormap
numColors = 256; % Typical colormap length

% Create a matrix of RGB values that transition from yellow to red
% The red channel stays at 1, the green channel transitions from 1 to 0, and the blue channel stays at 0
customColormap = [ones(numColors, 1) linspace(1, 0, numColors)' zeros(numColors, 1)];


GREEN = shaperead('../shapefile/BerlinGreen');                              % Reading shapefile of green spaces
WATER = shaperead('../shapefile/BerlinWater');                              % Reading shapefile of water spaces
SHAPE = shaperead(shapefile);                                               % Obejct to refer to the user-specified shapefile
clf                                                                         % Clear any pre-existing figure 
SHOW = variable;                                                            % Object to refer to outcome chosen by user
for i = 1:size(SHOW,1)                                                      % We change the sign of the outcome so that larger vales correspond to red in the chosen colour scheme
   SHAPE(i).SHOW = SHOW(i);
end
SHOWmin = (min(SHOW))*1;                                                   % We read the minimum value of the outcome
SHOWmax = (max(SHOW))*1;                                                   % We read the max value of the outcome
colorRange = makesymbolspec('Polygon', ...                                  % We define some feature of the maps to be generated, such as that it is a polygon,
    {'SHOW',[ SHOWmin SHOWmax], 'FaceColor', colormap((customColormap))}, ...       % shows the object SHOW, uses the min and max value for the colour range, and what colour scheme to use
    {'Default','EdgeColor','none'});                                        % We do not show polygon outlines

Figure2 = mapshow(SHAPE, 'Symbolspec', colorRange);                         % Now we generage our figure using the specified shapefile
    hold on;                                
    mapshow(WATER, 'FaceColor',[0.678 0.847 0.902], 'EdgeColor', 'none');   % Add water shape
    mapshow(GREEN, 'FaceColor',[0.133 0.545 0.133], 'EdgeColor', 'none');   % Add green shape
    hcb = colorbar;                                                         % Next three lines add color bar.
    set(hcb, 'Ylim', [ SHOWmin SHOWmax]);
    set(get(hcb, 'Title'),'String','Value');                  
    title(name);                                                            % We use the user-specified title
    name1 =  '/';                                                           % The next lines generate a pathname and file name
    name2 = folder;
    name3 = '/';
  	name4 = filename;
    name5 = '.png';
    filename = [name1 name2 name3 name4 name5];
saveas(Figure2, [pwd filename]);                                            % The map is being saved under the chosen name at the chosen destination
RESULT = 'DONE';                                                            % Confirm successful execution
end

% Code ends 