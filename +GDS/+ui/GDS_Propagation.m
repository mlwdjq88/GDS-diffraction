classdef GDS_Propagation < mic.Base
    
    
    properties (Constant)
        dWidth  = 870;
        dHeight =  600;
        
        % Axes tab IDS
        U8FARFIELD          = 1
        U8PUPIL             = 2
        U8MASK              = 3
        U8ZERNIKE           = 4;
    end
    
    properties
        
        cAppPath = fileparts(mfilename('fullpath'))
        
        % Graphical elements
        hFigure     % Main figure (not overwritable)
        
        % valuables
        dGDSData
        dGDSDataHead
        dGDSDataEnd
        dPolygon
        dMask
        dFarfield
        dPupil
        dPupilAmp
        dPhase
        dRMS
        dZrn
        dUx
        dUy
        x_um
        y_um
        yp_um
        xp_um
        % axis tab
        uitgAxesDisplay     % displays axes:field, pupil
        % Tab: field
        uitgFieldDisplay     % displays axes:field amp and pha
        haFieldAmp
        haFieldPha
        
        % Tab: pupil
        uitgPupilDisplay     % displays axes:pupil amp and pha and mask
        haPupilAmp
        haPupilPha
        haMask
        haZernike
        uitRMS
        
        % parameters
        hpPara
        uieLambda
        uieSampling
        uiePropDistance
        uieObjDistance
        uieOffsetAngle
        uieAzimuth
        uieRangeX
        uieRangeY
        uipPropMethod
        
        % control
        hpControl
        uieFilePath
        uibLoadFile
        uibPropagate
        
        % pupil phase analysis
        hpAnalysis
        uicbUnwrapPhase
        uicbRemoveSphere
        uicbRemoveTilt
        uicbRemoveDefocus
        uieNA
        uipMaskSelection
        uieNZernike
        uibDecomposition
        uibAnalyze
        
    end
    
    properties (SetAccess = private)
        
    end
    
    methods
        function this = GDS_Propagation()
            this.init()
        end
        
        
        
        function init(this)
            
            % axis tab
            this.uitgAxesDisplay = ...
                mic.ui.common.Tabgroup('ceTabNames', {'Far field', 'Pupil'});
            % Tab: field
            this.uitgFieldDisplay = ...
                mic.ui.common.Tabgroup('ceTabNames', {'Field amplitude', 'Field phase'});
            % Tab: pupil
            this.uitgPupilDisplay = ...
                mic.ui.common.Tabgroup('ceTabNames', {'Pupil amplitude', 'Pupil phase','Mask','Zernike'});
            this.uitRMS           = mic.ui.common.Text('cVal', 'RMS:');
            
            
            % parameters
            this.uieLambda       = mic.ui.common.Edit('cLabel', 'Wavelength(nm)', 'cType', 'd');
            this.uieSampling     = mic.ui.common.Edit('cLabel', 'Sampling(n¡Án)', 'cType', 'd');
            this.uiePropDistance = mic.ui.common.Edit('cLabel', 'Prop. distance(mm)', 'cType', 'd');
            this.uieObjDistance  = mic.ui.common.Edit('cLabel', 'Obj. distance(mm)', 'cType', 'd');
            this.uieOffsetAngle  = mic.ui.common.Edit('cLabel', 'Offset angle(¡ã)', 'cType', 'd');
            this.uieAzimuth      = mic.ui.common.Edit('cLabel', 'Azimuth(¡ã)', 'cType', 'd');
            this.uieRangeX         = mic.ui.common.Edit('cLabel', 'Field range X(um)', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieRangeY         = mic.ui.common.Edit('cLabel', 'Field range Y(um)', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uipPropMethod    = mic.ui.common.Popup('cLabel', 'Propagation method', 'ceOptions',...
                {'Polygon FFT','Fresnel integration'}, 'lShowLabel', true);
            this.uieLambda.set(13.5);
            this.uieSampling.set(32);
            this.uiePropDistance.set(1);
            this.uieObjDistance.set(0);
            this.uieOffsetAngle.set(0);
            this.uieAzimuth.set(0);
            this.uieRangeX.set('[-1, 1]');
            this.uieRangeY.set('[-1, 1]');
            this.uipPropMethod.setSelectedIndex(uint8(1));
            
            % control
            this.uieFilePath    = mic.ui.common.Edit('cLabel', 'GDS file path', 'cType', 'c','fhDirectCallback', @(src, evt)this.cb(src));
            this.uibLoadFile    = mic.ui.common.Button('cText', 'Load GDS', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibPropagate   = mic.ui.common.Button('cText', 'Propagate', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieFilePath.set('');
            
            % pupil phase analysis
            this.uicbUnwrapPhase   = mic.ui.common.Checkbox('cLabel', 'Unwrap phase');
            this.uicbRemoveSphere  = mic.ui.common.Checkbox('cLabel', 'Remove sphere');
            this.uicbRemoveTilt    = mic.ui.common.Checkbox('cLabel', 'Remove tilt');
            this.uicbRemoveDefocus = mic.ui.common.Checkbox('cLabel', 'Remove defocus');
            this.uieNA            = mic.ui.common.Edit('cLabel', 'NA', 'cType', 'd');
            this.uipMaskSelection = mic.ui.common.Popup('cLabel', 'Select mask', 'ceOptions',...
                {'No mask','Compute from pupil','Load mask'}, 'lShowLabel', true,'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieNZernike      = mic.ui.common.Edit('cLabel', 'Zernike #', 'cType', 'd');
            this.uibDecomposition = mic.ui.common.Button('cText', 'Decomposition', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibAnalyze = mic.ui.common.Button('cText', 'Analyze', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbUnwrapPhase.set(true);
            this.uicbRemoveSphere.set(false);
            this.uicbRemoveTilt.set(false);
            this.uicbRemoveDefocus.set(false);
            this.uieNA.set(0.0875);
            this.uieNZernike.set(24);
            this.uipMaskSelection.setSelectedIndex(uint8(1));
            
        end
        
        % Callback handler
        function cb(this, src,evt)
            switch src
                case {this.uieRangeX, this.uieRangeY}
                    this.validateCouplesEditBox(src, '[-1, 1]');
                    
                case this.uieFilePath
                    filename = this.uieFilePath.get();
                    if ~isempty(filename)
                        if ~strcmp(filename,'Please load GDS file first!')
                            this.dataLoading(filename);
                        end
                    else
                        this.uieFilePath.set('Please load GDS file first!');
                    end
                    
                case this.uibLoadFile
                    cDataDir = fullfile(this.cAppPath, '..','..','data', '*.gds');
                    [d, p] = uigetfile(cDataDir);
                    if isequal(d,0)||isequal(p,0)
                        return;
                    end
                    filename = [p d];
                    this.uieFilePath.set(filename);
                    this.dataLoading(filename);
                    
                case this.uibPropagate
                    this.propagate();
                    
                case this.uipMaskSelection
                    id =this.uipMaskSelection.getSelectedIndex();
                    switch id
                        case 1
                            this.dMask = [];
                        case 2
                            pupilAmp = abs(this.dPupil);
                            this.dMask = ones(size(this.dPupil));
                            this.dMask(pupilAmp<max(pupilAmp(:))/10) = 0;
                        case 3
                            cDataDir = fullfile(this.cAppPath, '..','..','data','mask', '*.mat');
                            [fn, pn] = uigetfile(cDataDir);
                            if isequal(fn,0)||isequal(pn,0)
                                return;
                            end
                            fileformat=fn(end-2:end);
                            filename= strcat(pn,fn);
                            switch fileformat
                                case 'mat'
                                    load(filename);
                            end
                            this.dMask = mask;
                    end
                    
                    % active tab
                    this.uitgAxesDisplay.selectTabByIndex(this.U8PUPIL);
                    this.uitgPupilDisplay.selectTabByIndex(this.U8MASK);
                    this.replot(this.U8MASK);
                    
                case this.uibDecomposition
                    this.decomposition();
                    % active tab
                    this.uitgAxesDisplay.selectTabByIndex(this.U8PUPIL);
                    this.uitgPupilDisplay.selectTabByIndex(this.U8ZERNIKE);
                    this.replot(this.U8ZERNIKE);
                    
                case this.uibAnalyze
                    this.analyze();
                    % active tab
                    this.uitgAxesDisplay.selectTabByIndex(this.U8PUPIL);
                    this.uitgPupilDisplay.selectTabByIndex(this.U8PUPIL);
                    this.replot(this.U8PUPIL);
                    
            end
        end
        
        function dataLoading(this, filename)
            % load data
            tic,
            outputFile=fopen(filename,'rb');
            if outputFile==-1
                fprintf('Opening file failed!\n');
                return;
            end
            offset=0;
            while 1
                fseek(outputFile, offset, 'bof');
                temp=fread(outputFile, [4,1],'uint8');
                if temp(3)==8&&temp(4)==0
                    break;
                else
                    offset=temp(1)*16+temp(2)+offset;
                end
            end
            fseek(outputFile, offset, 'bof');
            gds=fread(outputFile, 'int32','b');
            gds(end-1:end)=[];
            s1=find(gds==264192);% head of the polygon
            s2=find(gds==266496);% end of the polygon
            % remove error datas
            ds1 = diff(s1);
            ds2 = diff(s2);
            for i = 1:length(ds1)
                try
                    while ds1(i)>ds2(i)
                        ds2(i) = ds2(i)+ds2(i+1);
                        ds2(i+1) =[];
                        s2(i+1)  =[];
                    end
                    while ds1(i)<ds2(i)
                        ds1(i) = ds1(i)+ds1(i+1);
                        ds1(i+1) =[];
                        s1(i+1)  =[];
                    end
                catch
                    break;
                end
            end
            num=length(s1);
            fclose(outputFile);
            if num~=length(s2)
                fprintf('Data format is not correct!\n');
                return;
            end
            this.dGDSData = gds;
            this.dGDSDataHead = s1;
            this.dGDSDataEnd = s2;
            fprintf('Reading GDS file took %ds\n',round(toc));
        end
        
        % propagate function
        function propagate(this)
            % --  Exposure tool information
            s2 = this.dGDSDataEnd;
            s1 = this.dGDSDataHead;
            gds = this.dGDSData;
            if isempty(gds)
                fprintf('Please reload datafile!\n');
                return;
            end
            tic,
            num=length(s1);
            wavl = this.uieLambda.get();
            nrd = this.uieSampling.get();
            offsetAngle = this.uieOffsetAngle.get()*pi/180;
            azimuth = this.uieAzimuth.get()*pi/180;
            propdis_nm=this.uiePropDistance.get()*1e6;
            defocus_mm = this.uieObjDistance.get();
            propMethod = this.uipPropMethod.getSelectedIndex();
            if propMethod ==1
                sph = @(x,y)(2*pi/wavl*sign(defocus_mm)*sqrt(x.^2+y.^2+...
                    (defocus_mm*1e6+x*tan(offsetAngle)*cos(azimuth)+y*tan(offsetAngle)*sin(azimuth)).^2)+...% point source illumination
                    2*pi/wavl*sqrt(x.^2+y.^2+...
                    (propdis_nm-x*tan(offsetAngle)*cos(azimuth)-y*tan(offsetAngle)*sin(azimuth)).^2)); % compensate phase;
            else
                sph = @(x,y)(2*pi/wavl*sign(defocus_mm)*sqrt(x.^2+y.^2+...
                    (defocus_mm*1e6+x*tan(offsetAngle)*cos(azimuth)+y*tan(offsetAngle)*sin(azimuth)).^2));% point source illumination
            end
            rangeX = eval(this.uieRangeX.get());
            rangeY = eval(this.uieRangeY.get());
            xLeft= rangeX(1);
            xRight= rangeX(2);
            yLeft= rangeY(1);
            yRight= rangeY(2);
            
            xc=linspace(sin(atan(xLeft/propdis_nm*1000)),sin(atan(xRight/propdis_nm*1000)),2*nrd-1)/wavl;
            yc=linspace(sin(atan(yLeft/propdis_nm*1000)),sin(atan(yRight/propdis_nm*1000)),2*nrd-1)/wavl;
            x_um=propdis_nm*tan(asin(xc*wavl))/1000;
            y_um=propdis_nm*tan(asin(yc*wavl))/1000;
            [x_nm,y_nm] = meshgrid(x_um*1000,y_um*1000);
            % [x_nm,y_nm] = meshgrid(linspace(xLeft,xRight,2*nrd-1)*1000,linspace(yLeft,yRight,2*nrd-1)*1000);
            [nyux,nyuy]=meshgrid(xc,yc); % frequency coordinates
            this.dUx = nyux;
            this.dUy = nyuy;
            this.x_um = x_um;
            this.y_um = y_um;
            dpx_um = wavl*propdis_nm/(xRight-xLeft)*1e-6;
            dpy_um = wavl*propdis_nm/(yRight-yLeft)*1e-6;
            this.xp_um = dpx_um*[-nrd+1:nrd-1];
            this.yp_um = dpy_um*[-nrd+1:nrd-1];
            %% create polygon structs
            polyg(num,1)=struct('xy',[],'tx',[],'phase',[]);
            backg.int=0;
            backg.phase=0;
            for j=1:num
                if abs(cos(azimuth))==1
                    xs=gds((s1(j)+5):2:(s2(j)-3))/10*cos(offsetAngle); % from A to nm (tilt only works for x axis)
                    ys=gds((s1(j)+6):2:(s2(j)-3))/10;
                else
                    xs=gds((s1(j)+5):2:(s2(j)-3))/10; % from A to nm (tilt only works for x axis)
                    ys=gds((s1(j)+6):2:(s2(j)-3))/10*cos(offsetAngle);
                end
                polyg(j).xy =[xs';ys'];
                polyg(j).tx = 1;
                x0=mean(xs);
                y0=mean(ys);
                polyg(j).phase = sph(x0,y0);
                %     figure(2),plot3(x0,y0,sph(x0,y0),'.'),hold on;
            end
            
            % Forces all polygons to be defined in a CW orientation
            polyg = GDS.utils.cc_or_ccw(polyg) ;
            
            % -- Pupil intensity calculation
            % Computes the fourier transform of polygons in the plane specified by
            if propMethod ==1
                this.dFarfield = GDS.utils.polyProp(polyg, backg, nrd,nyux,nyuy) ;
            else
                this.dFarfield = GDS.utils.sphereProp(polyg,x_nm,y_nm,propdis_nm,wavl,offsetAngle,azimuth);
            end
            this.dPolygon = polyg;
            fprintf('Propagation took %ds\n',round(toc));
            % Make Field tab active:
            this.uitgAxesDisplay.selectTabByIndex(this.U8FARFIELD);
            this.replot(this.U8FARFIELD);
        end
        
        function analyze(this)
            diff_amp = this.dFarfield;
            norm_flag = 1 ; % 1: Normalize the intensity, 0: Unnormalize the intensity
            pupil = ifftshift(ifft2(ifftshift(diff_amp)));
            this.dPupil = pupil;
            if norm_flag>0
                pupil = pupil./max(abs(pupil(:))) ;
            else
                pupil = pupil.*length(pupil).^2 ;
            end
            pupilAmp = abs(pupil);
            pupilPha = atan2(imag(pupil),real(pupil));
            mask = this.dMask;
            if length(mask)~=length(pupil)||isempty(mask)
                mask = ones(length(pupil));
            end
            if this.uicbUnwrapPhase.get()
                pupilPha = GDS.utils.UnwrapPhaseBySortingReliabilityWithMask(pupilPha,mask*255);
            end
            pupilPha(mask==0) = NaN;
            pupilAmp(mask==0) = NaN;
            if this.uicbRemoveSphere.get()
                propdis_um=this.uiePropDistance.get()*1e3;
                NA = sin(atan((max(this.xp_um)-min(this.xp_um))/2/propdis_um));
                pupilPha = GDS.utils.removeSphere(pupilPha,NA,mask);
            end
            if this.uicbRemoveTilt.get()
                pupilPha = GDS.utils.DelTilt(pupilPha);
            end
            if this.uicbRemoveDefocus.get()
                pupilPha = GDS.utils.DelDefocus(pupilPha);
            end
            this.dPhase = pupilPha;
            this.dPupilAmp = pupilAmp;
            this.dRMS = std(pupilPha(mask==1))/2/pi;
        end
        
        function decomposition(this)
            pupilPha = this.dPhase;
            if isempty(pupilPha)
                this.analyze();
                pupilPha = this.dPhase;
            end
            NA = this.uieNA.get();
            propdis_um=this.uiePropDistance.get()*1e3;
            R_um = propdis_um*tan(asin(NA));
            [x,y] = meshgrid(this.xp_um,this.yp_um);
            [sx,sy] = size(pupilPha);
            order =this.uieNZernike.get();
            X=x/R_um;
            Y=y/R_um;
            basis = zeros(sx*sy, order + 1);
            mask = this.dMask;
            if isempty(mask)
                mask = ones(size(pupilPha));
            end
            basis(:,1) = mask(:);
            for k = 1:order
                [n, m] = GDS.utils.j2nm(k);
                zRTH = GDS.utils.ZgenNM(n,m);
                [TH, R] = cart2pol(X(:), Y(:));
                basis(:, k+1) = zRTH(R,TH).*mask(:);
            end
            this.dZrn =pinv(basis(mask==1,:))*pupilPha(mask==1);
            
        end
        
        
        % validates whether a char edit box evaluates to a Nx2 matrix,
        % colors accordingly.  Empty value is changed to []
        function [lOut, vals] = validateCouplesEditBox(~, src, cDefaultVal)
            lOut = true;
            vals = [];
            if isempty(src.get())
                src.styleDefault();
                src.set(cDefaultVal);
                return
            end
            try
                vals = eval(src.get());
                [~, sc] = size(vals);
                if (sc == 2 || sc == 0)
                    src.styleDefault();
                    lOut = true;
                else
                    src.styleBad();
                    lOut = false;
                end
            catch
                % can't read this edit box
                src.styleBad();
                lOut = false;
            end
        end
        
        % Main redraw function. Pass tab indices to refresh axes
        function replot(this, dTabIdx)
            
            switch dTabIdx
                
                case this.U8FARFIELD
                    diff_amp = this.dFarfield;
                    norm_flag = 1 ; % 1: Normalize the intensity, 0: Unnormalize the intensity
                    if norm_flag>0
                        diff_amp = diff_amp./max(abs(diff_amp(:))) ;
                    else
                        diff_amp = diff_amp./length(diff_amp).^2 ;
                    end
                    fieldAmp=abs(diff_amp);
                    fieldPha=atan2(imag(diff_amp),real(diff_amp));
                    
                    imagesc(this.haFieldAmp,this.x_um,this.y_um,fieldAmp);axis(this.haFieldAmp,'xy'); colorbar(this.haFieldAmp);
                    xlabel(this.haFieldAmp,'x/um'),ylabel(this.haFieldAmp,'y/um');
                    imagesc(this.haFieldPha,this.x_um,this.y_um,fieldPha/2/pi);axis(this.haFieldPha,'xy'); colorbar(this.haFieldPha);
                    xlabel(this.haFieldPha,'x/um'),ylabel(this.haFieldPha,'y/um');
                    
                case this.U8PUPIL
                    pupilAmp = this.dPupilAmp;
                    pupilPha = this.dPhase;
                    this.uitRMS.set(['RMS: ',num2str(round(this.dRMS*1000)/1000),' waves']);
                    imagesc(this.haPupilAmp,this.xp_um,this.yp_um,pupilAmp);axis(this.haPupilAmp,'xy'); colorbar(this.haPupilAmp);
                    xlabel(this.haPupilAmp,'x/um'),ylabel(this.haPupilAmp,'y/um');
                    imagesc(this.haPupilPha,this.xp_um,this.yp_um,pupilPha/2/pi);axis(this.haPupilPha,'xy'); colorbar(this.haPupilPha);
                    xlabel(this.haPupilPha,'x/um'),ylabel(this.haPupilPha,'y/um');
                    
                case this.U8MASK
                    imagesc(this.haMask,this.dMask);axis(this.haMask,'xy','off','equal','tight');
                    
                case this.U8ZERNIKE
                    bar(this.haZernike,0:(this.uieNZernike.get()),this.dZrn);
                    xlabel(this.haZernike,'Zernike terms'),ylabel(this.haZernike,'Coefficients / Waves');
            end
            
        end
        
        
        
        
        function build(this, hFigure, dOffsetX, dOffsetY)
            if nargin <3
                dOffsetX = 0;
                dOffsetY = 0;
            elseif nargin == 3
                dOffsetY = dOffsetX;
                dOffsetX =  hFigure;            
            end
            
            % build the main window
            if nargin == 2||nargin == 4
                this.hFigure = hFigure;
            else
                this.hFigure = figure(...
                    'name', 'GDS propagation GUI v1.200701',...
                    'Units', 'pixels',...
                    'Position', [5 - dOffsetX, 5 - dOffsetY,  this.dWidth, this.dHeight],...
                    'handlevisibility','off',... %out of reach gcf
                    'numberTitle','off',...
                    'Toolbar','none',...
                    'Menubar','none');
            end
            
            
            % Build all containers first:
            drawnow
            
            % Axes
            dTgPx = 20;
            dTgPy = 20;
            this.uitgAxesDisplay.build(this.hFigure, dTgPx, dTgPy, 550, 550);
            
            % Axes:Far field
            uitField = this.uitgAxesDisplay.getTabByName('Far field');
            this.uitgFieldDisplay.build(uitField, dTgPx, dTgPy+20, 510, 500);
            uitFieldAmp = this.uitgFieldDisplay.getTabByName('Field amplitude');
            
            this.haFieldAmp = axes('Parent', uitFieldAmp, ...
                'Units', 'pixels', ...
                'Position', [50, 60, 430, 360], ...
                'XTick', [], 'YTick', []);
            uitFieldPha = this.uitgFieldDisplay.getTabByName('Field phase');
            
            this.haFieldPha = axes('Parent', uitFieldPha, ...
                'Units', 'pixels', ...
                'Position', [50, 60, 430, 360], ...
                'XTick', [], 'YTick', []);
            
            % Axes:pupil
            uitPupil = this.uitgAxesDisplay.getTabByName('Pupil');
            
            this.uitgPupilDisplay.build(uitPupil, dTgPx, dTgPy+20, 510, 500);
            uitPupilAmp = this.uitgPupilDisplay.getTabByName('Pupil amplitude');
            this.haPupilAmp = axes('Parent', uitPupilAmp, ...
                'Units', 'pixels', ...
                'Position', [50, 60, 430, 360], ...
                'XTick', [], 'YTick', []);
            uitPupilPha = this.uitgPupilDisplay.getTabByName('Pupil phase');
            this.haPupilPha = axes('Parent', uitPupilPha, ...
                'Units', 'pixels', ...
                'Position', [50, 60, 430, 360], ...
                'XTick', [], 'YTick', []);
            uitMask = this.uitgPupilDisplay.getTabByName('Mask');
            this.haMask = axes('Parent', uitMask, ...
                'Units', 'pixels', ...
                'Position', [40, 60, 430, 360], ...
                'XTick', [], 'YTick', []);
            uitZernike = this.uitgPupilDisplay.getTabByName('Zernike');
            this.haZernike = axes('Parent', uitZernike, ...
                'Units', 'pixels', ...
                'Position', [55, 60, 430, 360], ...
                'XTick', [], 'YTick', []);
            
            this.uitRMS.build           (uitPupilPha, 80, 40, 200, 20);
            this.uitRMS.setFontSize(14);
            
            % parameters
            this.hpPara = uipanel(...
                'Parent', this.hFigure,...
                'Units', 'pixels',...
                'Title', 'Parameters',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [585 350 260 230] ...
                );
            this.uieLambda.build (this.hpPara,20,20,100,20);
            this.uieSampling.build (this.hpPara,140,20,100,20);
            this.uiePropDistance.build (this.hpPara,20,60,100,20);
            this.uieObjDistance.build (this.hpPara,140,60,100,20);
            this.uieOffsetAngle.build (this.hpPara,20,100,100,20);
            this.uieAzimuth.build (this.hpPara,140,100,100,20);
            this.uieRangeX.build (this.hpPara,20,140,100,20);
            this.uieRangeY.build (this.hpPara,140,140,100,20);
            this.uipPropMethod.build (this.hpPara,20,180,100,20);
            
            % control
            this.hpControl = uipanel(...
                'Parent', this.hFigure,...
                'Units', 'pixels',...
                'Title', 'Control',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [585 240 260 100] ...
                );
            this.uieFilePath.build (this.hpControl,20,20,220,20);
            this.uibLoadFile.build (this.hpControl,20,70,100,20);
            this.uibPropagate.build (this.hpControl,140,70,100,20);
            
            % pupil phase analysis
            this.hpAnalysis = uipanel(...
                'Parent', this.hFigure,...
                'Units', 'pixels',...
                'Title', 'Analysis',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [585 30 260 200] ...
                );
            this.uipMaskSelection.build (this.hpAnalysis,20,20,220,20);
            this.uicbUnwrapPhase.build (this.hpAnalysis,20,70,100,20);
            this.uicbRemoveSphere.build (this.hpAnalysis,140,70,100,20);
            this.uicbRemoveTilt.build (this.hpAnalysis,20,100,100,20);
            this.uicbRemoveDefocus.build (this.hpAnalysis,140,100,110,20);
            this.uieNA.build (this.hpAnalysis,140,125,100,20);
            this.uieNZernike.build (this.hpAnalysis,20,125,100,20);
            this.uibDecomposition.build (this.hpAnalysis,20,170,100,20);
            this.uibAnalyze.build (this.hpAnalysis,140,170,100,20);
            drawnow;
        end
        
    end
    
end

