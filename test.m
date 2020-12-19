%% test script
if ~exist('gds')||~ishandle(gds.hFigure)
    launch_GDS;
end

%% parameter setting
gds.uieLambda.set(13.56);
gds.uieSampling.set(64);
gds.uiePropDistance.set(0.4);
gds.uieObjDistance.set(43.4);
gds.uieOffsetAngle.set(6);
gds.uieAzimuth.set(180);
gds.uieRangeX.set('[-4, 4]');
gds.uieRangeY.set('[-4, 4]');
gds.uipPropMethod.setSelectedIndex(uint8(1));

% file path, select gds file here
gds.uieFilePath.set('D:\OneDrive\GDS-diffraction\data\SERM_400_correctFlip_wrv_v3.gds');
drawnow;

% propagate
gds.cb(gds.uibPropagate);

%% analysis para
gds.uicbUnwrapPhase.set(true);
gds.uicbRemoveSphere.set(false);
gds.uicbRemoveTilt.set(false);
gds.uicbRemoveDefocus.set(false);
gds.uieNA.set(0.0875);
gds.uieNZernike.set(24);
gds.uipMaskSelection.setSelectedIndex(uint8(2));
gds.cb(gds.uipMaskSelection);

%% load mask
% gds.uipMaskSelection.setSelectedIndex(uint8(3));
% gds.cb(gds.uipMaskSelection);
%% pupil analysis
gds.cb(gds.uibAnalyze);

%% Zernike decomposition
gds.cb(gds.uibDecomposition);



