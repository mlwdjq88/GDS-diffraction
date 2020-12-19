addpath('../mpm');

% mic library
mpm clearpath
mpm addpath

% Add this root directory to path
[cDirThis, cName, cExt] = fileparts(mfilename('fullpath'));
addpath(genpath(cDirThis));



% Build LSI UI
hFigure = figure(...
    'name', 'GDS propagation GUI v1.200701',...
    'Units', 'pixels',...
    'Position', [500, 300,  900, 600],...
    'handlevisibility','off',... %out of reach gcf
    'numberTitle','off',...
    'Toolbar','none',...
    'Menubar','none');

gds = GDS.ui.GDS_Propagation;
gds.build(hFigure,-500,-300);
