function varargout = GDSDiff(varargin)
% GDSDIFF MATLAB code for GDSDiff.fig
%      GDSDIFF, by itself, creates a new GDSDIFF or raises the existing
%      singleton*.
%
%      H = GDSDIFF returns the handle to a new GDSDIFF or the handle to
%      the existing singleton*.
%
%      GDSDIFF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GDSDIFF.M with the given input arguments.
%
%      GDSDIFF('Property','Value',...) creates a new GDSDIFF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GDSDiff_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GDSDiff_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GDSDiff

% Last Modified by GUIDE v2.5 29-Jun-2020 16:19:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GDSDiff_OpeningFcn, ...
                   'gui_OutputFcn',  @GDSDiff_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GDSDiff is made visible.
function GDSDiff_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GDSDiff (see VARARGIN)

% Choose default command line output for GDSDiff
handles.output = hObject;
try
    mpm addpath
catch
    addpath('mpm-packages\ryan_toolbox');
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GDSDiff wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GDSDiff_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% parpool;
% Get default command line output from handles structure
varargout{1} = handles.output;



function GDSPN_Callback(hObject, eventdata, handles)
% hObject    handle to GDSPN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GDSPN as text
%        str2double(get(hObject,'String')) returns contents of GDSPN as a double


% --- Executes during object creation, after setting all properties.
function GDSPN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GDSPN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadGDS.
function loadGDS_Callback(hObject, eventdata, handles)
% hObject    handle to loadGDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global old_pr
if old_pr==0
    old_pr=[];
end
[fn,pn]=uigetfile({'*.gds','GDS file';'*.*','All files'},'Load GDS file',old_pr);
if isequal(fn,0)||isequal(pn,0)
    return;
end
old_pr=pn;
GDSPN=strcat(pn,fn);
set(handles.GDSPN,'String',GDSPN);



function fieldSampling_Callback(hObject, eventdata, handles)
% hObject    handle to fieldSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldSampling as text
%        str2double(get(hObject,'String')) returns contents of fieldSampling as a double


% --- Executes during object creation, after setting all properties.
function fieldSampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function distance_Callback(hObject, eventdata, handles)
% hObject    handle to distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance as text
%        str2double(get(hObject,'String')) returns contents of distance as a double


% --- Executes during object creation, after setting all properties.
function distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xLeft_Callback(hObject, eventdata, handles)
% hObject    handle to xLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xLeft as text
%        str2double(get(hObject,'String')) returns contents of xLeft as a double


% --- Executes during object creation, after setting all properties.
function xLeft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xRight_Callback(hObject, eventdata, handles)
% hObject    handle to xRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xRight as text
%        str2double(get(hObject,'String')) returns contents of xRight as a double


% --- Executes during object creation, after setting all properties.
function xRight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yLeft_Callback(hObject, eventdata, handles)
% hObject    handle to yLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yLeft as text
%        str2double(get(hObject,'String')) returns contents of yLeft as a double


% --- Executes during object creation, after setting all properties.
function yLeft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yRight_Callback(hObject, eventdata, handles)
% hObject    handle to yRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yRight as text
%        str2double(get(hObject,'String')) returns contents of yRight as a double


% --- Executes during object creation, after setting all properties.
function yRight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function waveLength_Callback(hObject, eventdata, handles)
% hObject    handle to waveLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of waveLength as text
%        str2double(get(hObject,'String')) returns contents of waveLength as a double


% --- Executes during object creation, after setting all properties.
function waveLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to waveLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%read GDS
if get(hObject,'Value')==0
    set(hObject,'String','Run');
    return;
else
    set(hObject,'String','Stop');
end
drawnow;
filename=get(handles.GDSPN,'String');
outputFile=fopen(filename,'rb');
if outputFile==-1
    return;
end
% fseek(outputFile, 0, 'eof');
% uint8Length=ftell(outputFile);
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
    set(handles.run,'String','run');
    return;
end
tic;
% --  Exposure tool information
wavl = str2double(get(handles.waveLength,'String'));
nrd = str2double(get(handles.fieldSampling,'String'));
offsetAngle = str2double(get(handles.uieAngle,'String'))*pi/180;
azimuth = str2double(get(handles.uieAzimuth,'String'))*pi/180;
propdis_nm=str2double(get(handles.distance,'String'))*1e6;
defocus_mm = str2double(get(handles.uieDefocus,'String'));
propMethod = get(handles.uipPropMethod,'Value');
if propMethod ==1
    sph = @(x,y)(2*pi/wavl*sign(defocus_mm)*sqrt(x.^2+y.^2+...
        (defocus_mm*1e6+x*tan(offsetAngle)*cos(azimuth)+y*tan(offsetAngle)*sin(azimuth)).^2)+...% point source illumination
        2*pi/wavl*sqrt(x.^2+y.^2+...
        (propdis_nm-x*tan(offsetAngle)*cos(azimuth)-y*tan(offsetAngle)*sin(azimuth)).^2)); % compensate phase;
else
    sph = @(x,y)(2*pi/wavl*sign(defocus_mm)*sqrt(x.^2+y.^2+...
        (defocus_mm*1e6+x*tan(offsetAngle)*cos(azimuth)+y*tan(offsetAngle)*sin(azimuth)).^2));% point source illumination
end
xLeft= str2double(get(handles.xLeft,'String'));
xRight= str2double(get(handles.xRight,'String'));
yLeft= str2double(get(handles.yLeft,'String'));
yRight= str2double(get(handles.yRight,'String'));
xc=linspace(sin(atan(xLeft/propdis_nm*1000)),sin(atan(xRight/propdis_nm*1000)),2*nrd-1)/wavl;
yc=linspace(sin(atan(yLeft/propdis_nm*1000)),sin(atan(yRight/propdis_nm*1000)),2*nrd-1)/wavl;
x_um=propdis_nm*tan(asin(xc*wavl))/1000;
y_um=propdis_nm*tan(asin(yc*wavl))/1000;
[x_nm,y_nm] = meshgrid(x_um*1000,y_um*1000);
% [x_nm,y_nm] = meshgrid(linspace(xLeft,xRight,2*nrd-1)*1000,linspace(yLeft,yRight,2*nrd-1)*1000);
[nyux,nyuy]=meshgrid(xc,yc); % frequency coordinates
norm_flag = 1 ; % 1: Normalize the intensity, 0: Unnormalize the intensity
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
polyg = cc_or_ccw(polyg) ;

% -- Pupil intensity calculation
% Computes the fourier transform of polygons in the plane specified by
if propMethod ==1
    diff_amp = polyProp(polyg, backg, nrd,nyux,nyuy,handles) ;
else
    diff_amp = sphereProp(polyg,x_nm,y_nm,propdis_nm,wavl,offsetAngle,azimuth,handles);
end

pupil = ifftshift(ifft2(ifftshift(diff_amp)));
if norm_flag>0
    diff_amp = diff_amp./max(abs(diff_amp(:))) ;
    pupil = pupil./max(abs(pupil(:))) ;
else
    diff_amp = diff_amp./(2.*nrd-1).^2 ;
    pupil = pupil.*(2.*nrd-1).^2 ;
end
fieldAmp=abs(diff_amp);
fieldPha=atan2(imag(diff_amp),real(diff_amp));
pupilAmp = abs(pupil);
pupilPha = atan2(imag(pupil),real(pupil));

dpx_um = wavl*propdis_nm/(xRight-xLeft)*1e-6;
dpy_um = wavl*propdis_nm/(yRight-yLeft)*1e-6;
xp_um = dpx_um*[-nrd+1:nrd-1];
yp_um = dpy_um*[-nrd+1:nrd-1];
imagesc(handles.outputField,x_um,y_um,fieldAmp); colorbar(handles.outputField);
xlabel(handles.outputField,'x/um'),ylabel(handles.outputField,'y/um');title(handles.outputField,'Field amplitude');
imagesc(handles.outputPhase,x_um,y_um,fieldPha/2/pi); colorbar(handles.outputPhase);
xlabel(handles.outputPhase,'x/um'),ylabel(handles.outputPhase,'y/um');title(handles.outputPhase,'Field phase');
imagesc(handles.pupilAmp,xp_um,yp_um,pupilAmp); colorbar(handles.pupilAmp);
xlabel(handles.pupilAmp,'x/um'),ylabel(handles.pupilAmp,'y/um');title(handles.pupilAmp,'Pupil amplitude');
imagesc(handles.pupilPhase,xp_um,yp_um,pupilPha/2/pi); colorbar(handles.pupilPhase);
xlabel(handles.pupilPhase,'x/um'),ylabel(handles.pupilPhase,'y/um');title(handles.pupilPhase,'Pupil phase');
set(hObject,'String','Run','Value',0);
setappdata(gcf,'field',diff_amp);
setappdata(gcf,'pupil',pupil);
setappdata(gcf,'pupilPha',pupilPha);
setappdata(gcf,'x_um',x_um);
setappdata(gcf,'y_um',y_um);
setappdata(gcf,'xp_um',xp_um);
setappdata(gcf,'yp_um',yp_um);
fprintf('Propagation took %0.1fs\n',toc);

function p_amp = sphereProp(polyg,x_nm,y_nm,z_nm,lambda_nm,offsetAngle,azimuth,handles) 
p_amp =zeros(size(x_nm));
num_poly = length(polyg);
runtimes=datenum(datestr(now,31));
parfor i_poly = 1:num_poly
     xy0 = mean(polyg(i_poly).xy,2);
     x = xy0(1);
     y = xy0(2);
     z =x*tan(offsetAngle)*cos(azimuth)+y*tan(offsetAngle)*sin(azimuth);
     R_nm = sqrt((x-x_nm).^2+(y-y_nm).^2 +(z-z_nm).^2);
     p_amp = p_amp + polyg(i_poly).tx.*exp(1i*polyg(i_poly).phase+1i*2*pi/lambda_nm*R_nm)./R_nm.^2*(z_nm-z);
end
rt=datestr(datenum(datestr(now,31))-runtimes,13);
set(handles.RT,'String',rt);


function p_amp = polyProp(polyg,backg,nrd,nyux,nyuy,handles) 
% Function M-file polyFFT.m
% *************************************************************************
% Modified for the University of California: 5/9/2010
% Last modified date: 10/2/2014
%
% - Use this code at your own risk.
% - Freely used inside the University of California.
% - Do not use commercially.
% - If you make a significant improvement, please let me know to 
%   kenji.yamazoe@gmail.com.
%
% *************************************************************************
% ****** Copyright (C) Kenji Yamazoe, 2004-2014. All Right Reserved. ******
% *************************************************************************

% II=sqrt(-1) ;

% region = (wavl./(2.*NA).*(nrd-1)) ;

% Generate coordinates of observation plane
% nyu = (-(nrd-1):(nrd-1)).*(2./(nrd-1)).*NA./wavl ;
% [nyux,nyuy] = meshgrid(nyu) ;


% generate FFT of mask aperture
% p_amp = sqrt(backg.int).*exp(II.*backg.phase) ...
%       .* region.^2.*feval(@sinc,region.*nyux).*feval(@sinc,region.*nyuy) ;
runtimes=datenum(datestr(now,31));
p_amp =zeros(size(nyux));

% The number of the polygons
num_poly = length(polyg) ;
oldelapsedTime=100;
xunit_vec = [1,0] ;
parfor i_poly = 1:num_poly
%     drawnow;
    stopsign=get(handles.run,'value');
    if stopsign==0
%          break;
    end
     timerVal=tic;
    dp_x = [polyg(i_poly).xy(1,1:end),polyg(i_poly).xy(1,1)] ;
    dp_y = [polyg(i_poly).xy(2,1:end),polyg(i_poly).xy(2,1)] ;
    
%     max_x = max(dp_x) ;
%     min_x = min(dp_x) ;
%     max_y = max(dp_y) ;
    min_y = min(dp_y) ;
    
    for i=1:(length(dp_x)-1)
        
        vec1 = [dp_x(i+1)-dp_x(i),dp_y(i+1)-dp_y(i)] ;
        direc = dot(vec1,xunit_vec)./norm(vec1) ;
        
        % Vertical line - no action
        % complete on 4/4/09
        if abs(direc)<eps
%              fprintf('Y-direction.\n') ;
            continue
        end
        
        % Holizontal line
        % complete on 4/4/09
        if abs(abs(direc)-1)<eps
%              fprintf('X-direction.\n') ;
            % If starting coordinate is > minimum y
            if abs(dp_y(i)-min_y)>eps
                p_amp = p_amp ...
                      + feval(@FFTrect,dp_x(i+1),dp_x(i),dp_y(i+1),min_y,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy) ;
            end
            continue
        end
        % diagonal line - left direction
        if direc<0
%              fprintf('Left direction: %5.3f. \n',direc) ;
            
            up_flag = 0 ;
            min_flag = 0 ;
            
            % If vector is going up
            if dp_y(i+1)>dp_y(i)
                up_flag = 1 ;
            end
            
            % If dp_y contains the minimum y value
            if abs(min([dp_y(i+1),dp_y(i)])-min_y)<eps
                min_flag = 1;
            end
            
            d_flag = up_flag.*2 + min_flag ;
            
            switch d_flag
                case 0
                    % left-down without minimun
                    tx = zeros(1,3) ;
                    ty = zeros(1,3) ;
                    tx(1) = dp_x(i+1) ;
                    tx(2) = dp_x(i) ;
                    tx(3) = dp_x(i) ;
                    ty(1) = dp_y(i+1) ;
                    ty(2) = dp_y(i+1) ;
                    ty(3) = dp_y(i) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTtri,tx,ty,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy,nrd) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTrect,dp_x(i+1),dp_x(i),dp_y(i+1),min_y,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy) ;

                case 1
                    % left-down with minimun
                    tx = zeros(1,3) ;
                    ty = zeros(1,3) ;
                    tx(1) = dp_x(i+1) ;
                    tx(2) = dp_x(i) ;
                    tx(3) = dp_x(i) ;
                    ty(1) = dp_y(i+1) ;
                    ty(2) = dp_y(i+1) ;
                    ty(3) = dp_y(i) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTtri,tx,ty,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy,nrd) ;
                    
                case 2
                    % left-up without minimun
                    tx = zeros(1,3) ;
                    ty = zeros(1,3) ;
                    tx(1) = dp_x(i+1) ;
                    tx(2) = dp_x(i) ;
                    tx(3) = dp_x(i+1) ;
                    ty(1) = dp_y(i) ;
                    ty(2) = dp_y(i) ;
                    ty(3) = dp_y(i+1) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTtri,tx,ty,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy,nrd) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTrect,dp_x(i+1),dp_x(i),dp_y(i),min_y,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy) ;
                      
                case 3
                    % left-up with minimun
                    tx = zeros(1,3) ;
                    ty = zeros(1,3) ;
                    tx(1) = dp_x(i+1) ;
                    tx(2) = dp_x(i) ;
                    tx(3) = dp_x(i+1) ;
                    ty(1) = dp_y(i) ;
                    ty(2) = dp_y(i) ;
                    ty(3) = dp_y(i+1) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTtri,tx,ty,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy,nrd) ;
                      
            end

            continue
        end
        
        % diagonal line - right direction
        if direc>0
%             fprintf('Right direction: %5.3f. \n',direc) ;
            
            up_flag = 0 ;
            min_flag = 0 ;
            
            if dp_y(i+1)>dp_y(i)
                up_flag = 1 ;
            end

            if abs(min([dp_y(i+1),dp_y(i)])-min_y)<eps
                min_flag = 1;
            end
            
            d_flag = up_flag.*2 + min_flag ;
            
            switch d_flag
                case 0
                    % right-down without minimun
                    tx = zeros(1,3) ;
                    ty = zeros(1,3) ;
                    tx(1) = dp_x(i) ;
                    tx(2) = dp_x(i+1) ;
                    tx(3) = dp_x(i) ;
                    ty(1) = dp_y(i) ;
                    ty(2) = dp_y(i+1) ;
                    ty(3) = dp_y(i+1) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTtri,tx,ty,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy,nrd) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTrect,dp_x(i+1),dp_x(i),dp_y(i+1),min_y,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy) ;

                case 1
                    % right-down with minimun
                    tx = zeros(1,3) ;
                    ty = zeros(1,3) ;
                    tx(1) = dp_x(i) ;
                    tx(2) = dp_x(i+1) ;
                    tx(3) = dp_x(i) ;
                    ty(1) = dp_y(i) ;
                    ty(2) = dp_y(i+1) ;
                    ty(3) = dp_y(i+1) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTtri,tx,ty,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy,nrd) ;
                    
                case 2
                    % right-up without minimun
                    tx = zeros(1,3) ;
                    ty = zeros(1,3) ;
                    tx(1) = dp_x(i) ;
                    tx(2) = dp_x(i+1) ;
                    tx(3) = dp_x(i+1) ;
                    ty(1) = dp_y(i) ;
                    ty(2) = dp_y(i+1) ;
                    ty(3) = dp_y(i) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTtri,tx,ty,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy,nrd) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTrect,dp_x(i+1),dp_x(i),dp_y(i),min_y,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy) ;
                      
                case 3
                    % right-up with minimun
                    tx = zeros(1,3) ;
                    ty = zeros(1,3) ;
                    tx(1) = dp_x(i) ;
                    tx(2) = dp_x(i+1) ;
                    tx(3) = dp_x(i+1) ;
                    ty(1) = dp_y(i) ;
                    ty(2) = dp_y(i+1) ;
                    ty(3) = dp_y(i) ;
                    
                    p_amp = p_amp ...
                          + feval(@FFTtri,tx,ty,polyg(i_poly).tx,polyg(i_poly).phase,backg.int,backg.phase,nyux,nyuy,nrd) ;
                    
            end

            continue

        end
        
%         fprintf('Nothing happened.\n') ;
        
    end
    elapsedTime= toc(timerVal);
%     if oldelapsedTime>elapsedTime
%         oldelapsedTime=elapsedTime;
%     else
%         elapsedTime=oldelapsedTime;
%     end
    rt=datestr(datenum(datestr(now,31))-runtimes,13);
    resttime=datestr(elapsedTime*(num_poly-i_poly)/24/3600+datenum('00:00:00','HH:MM:SS'),13);
    finishtime=datestr(datenum(datestr(now,31),'yyyy-mm-dd HH:MM:SS')+datenum(resttime)-datenum('00:00:00','HH:MM:SS'),31);
    if mod(i_poly,10)==0
        set(handles.RT,'String',rt);
        set(handles.ERT,'String',resttime);
        set(handles.EFT,'String',finishtime);
        drawnow;
    end
end
rt=datestr(datenum(datestr(now,31))-runtimes,13);
set(handles.RT,'String',rt);


% draw_t = linspace(0,2.*pi,101) ; 
% figure ; hold on ; 
% imagesc(nyu.*(wavl./NA),nyu.*(wavl./NA),real(p_amp)) ; axis equal ; colorbar ; plot(cos(draw_t),sin(draw_t),'LineWidth',2,'Color','w') ; axis equal tight ; axis([-2,2,-2,2]) ;
% set(gca,'TickDir','out','XTick',[-2,-1,0,1,2],'XTickLabel',[-2,-1,0,1,2],'YTick',[-2,-1,0,1,2],'YTickLabel',[-2,-1,0,1,2]) ;
% set(gca,'FontSize',16,'FontWeight','bold','FontName','Tahoma') ;
% xlabel('f','FontSize',20,'FontWeight','bold') ; ylabel('g','FontSize',20,'FontWeight','bold') ; 
% hold off
% 
% figure ; hold on ; 
% imagesc(nyu.*(wavl./NA),nyu.*(wavl./NA),imag(p_amp)) ; axis equal ; colorbar ; plot(cos(draw_t),sin(draw_t),'LineWidth',2,'Color','w') ; axis equal tight ; axis([-2,2,-2,2]) ;
% set(gca,'TickDir','out','XTick',[-2,-1,0,1,2],'XTickLabel',[-2,-1,0,1,2],'YTick',[-2,-1,0,1,2],'YTickLabel',[-2,-1,0,1,2]) ;
% set(gca,'FontSize',16,'FontWeight','bold','FontName','Tahoma') ;
% xlabel('f','FontSize',20,'FontWeight','bold') ; ylabel('g','FontSize',20,'FontWeight','bold') ; 
% hold off


% figure ; hold on ; 
% imagesc(propdis.*NA*nyu.*(wavl./NA),propdis.*NA*nyu.*(wavl./NA),abs(p_amp).^2) ; axis equal ; colorbar ;% plot(cos(draw_t),sin(draw_t),'LineWidth',2,'Color','w') ; %axis([-2,2,-2,2]) ;
% % set(gca,'TickDir','out','XTick',[-2,-1,0,1,2],'XTickLabel',[-2,-1,0,1,2],'YTick',[-2,-1,0,1,2],'YTickLabel',[-2,-1,0,1,2]) ;
%  axis equal tight ;
% set(gca,'FontSize',16,'FontWeight','bold','FontName','Tahoma') ;
% xlabel('X/um','FontSize',20,'FontWeight','bold') ; ylabel('Y/um','FontSize',20,'FontWeight','bold') ; 
% hold off

warning on MATLAB:divideByZero

% ------------------------------------------------------
% FFT of a rectangle
% ------------------------------------------------------
function y=FFTrect(x2,x1,y2,y1,tra_int,phase,tra_back,pha_back,nyux,nyuy)

II=sqrt(-1) ;

x = (x2+x1)./2 ;
y = (y2+y1)./2 ;
wx = x2 - x1 ;
wy = y2 - y1 ;

y = wx.*wy.*feval(@sinc,wx.*nyux).*feval(@sinc,wy.*nyuy).*exp(-II.*2.*pi.*(x.*nyux+y.*nyuy))...
              .*(sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)) ;

% ------------------------------------------------------
% FFT of a triangle
% ------------------------------------------------------
function y=FFTtri(tx,ty,tra_int,phase,tra_back,pha_back,nyux,nyuy,nrd)

II=sqrt(-1) ;

rot = 5.0e-8 ;
nyuxt = cos(rot).*nyux+sin(rot).*nyuy ;
nyuyt = -sin(rot).*nyux+cos(rot).*nyuy ;

% nyuxt = nyux ;
% nyuyt = nyuy ;
% 
% ty(1) = ty(1) - eps ;
% tx(2) = tx(2) - eps ;
% tx(3) = tx(3) + eps ;
% ty(3) = ty(3) + eps ;

% Formula
G_dom = 4.* pi.^2 ...
      .* (nyuxt .* (tx(1)-tx(2)) + nyuyt .* (ty(1)-ty(2))) ...
      .* (nyuxt .* (tx(3)-tx(1)) + nyuyt .* (ty(3)-ty(1))) ...
      .* (nyuxt .* (tx(3)-tx(2)) + nyuyt .* (ty(3)-ty(2))) ;
G_num = (tx(1).*(ty(3)-ty(2)) + tx(2).*(ty(1)-ty(3)) + tx(3).*(ty(2)-ty(1))) ...
      .*(exp(-II.*2.*pi.*(nyuxt.*tx(1) + nyuyt.*ty(1))) .* (nyuxt.*(tx(3)-tx(2)) + nyuyt.*(ty(3)-ty(2))) ...
       + exp(-II.*2.*pi.*(nyuxt.*tx(2) + nyuyt.*ty(2))) .* (nyuxt.*(tx(1)-tx(3)) + nyuyt.*(ty(1)-ty(3))) ...
       + exp(-II.*2.*pi.*(nyuxt.*tx(3) + nyuyt.*ty(3))) .* (nyuxt.*(tx(2)-tx(1)) + nyuyt.*(ty(2)-ty(1)))) ;

GT = (sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)).*G_num./G_dom ;
    
% formula in g = 0
if (abs(tx(2)-tx(3)))<3.*eps
    % this works well, 4/9/2009
    GT(nrd,:) = (tx(1).*(ty(3)-ty(2)) + tx(2).*(ty(1)-ty(3)) + tx(3).*(ty(2)-ty(1))) ...
              .*(II.*2.*pi.*nyux(1,:) .* exp(-II.*2.*pi.*nyux(1,:).*tx(3)).*(tx(3)-tx(1)) ...
               + exp(-II.*2.*pi.*nyux(1,:).*tx(3)) - exp(-II.*2.*pi.*nyux(1,:).*tx(1))) ...
              .*(sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)) ...
              ./ (4.*pi.^2.*nyux(1,:).^2 .* (tx(1)-tx(2)) .* (tx(1)-tx(3))) ;
end

if (abs(tx(1)-tx(3)))<3.*eps
    % this works well, 4/9/2009
    GT(nrd,:) = (tx(1).*(ty(3)-ty(2)) + tx(2).*(ty(1)-ty(3)) + tx(3).*(ty(2)-ty(1))) ...
              .*(II.*2.*pi.*nyux(1,:) .* exp(-II.*2.*pi.*nyux(1,:).*tx(3)).*(tx(2)-tx(3)) ...
               + exp(-II.*2.*pi.*nyux(1,:).*tx(2)) - exp(-II.*2.*pi.*nyux(1,:).*tx(3))) ...
              .*(sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)) ...
              ./ (4.*pi.^2.*nyux(1,:).^2 .* (tx(1)-tx(2)) .* (tx(2)-tx(3))) ;
end

if (abs(tx(1)-tx(2)))<3.*eps
    % I'm not sure if it works.
    GT(nrd,:) = (tx(1).*(ty(3)-ty(2)) + tx(2).*(ty(1)-ty(3)) + tx(3).*(ty(2)-ty(1))) ...
              .*(II.*2.*pi.*nyux(1,:) .* exp(-II.*2.*pi.*nyux(1,:).*tx(1)).*(tx(1)-tx(3)) ...
               + exp(-II.*2.*pi.*nyux(1,:).*tx(1)) - exp(-II.*2.*pi.*nyux(1,:).*tx(3))) ...
              .*(sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)) ...
              ./ (4.*pi.^2.*nyux(1,:).^2 .* (tx(2)-tx(3)) .* (tx(1)-tx(3))) ;
end

% formula in f = 0
if (abs(ty(2)-ty(3)))<3.*eps
    % this works well, 4/10/2009
    GT(:,nrd) = (tx(2).*(ty(1)-ty(3)) + tx(3).*(ty(2)-ty(1))) ...
              .*(II.*2.*pi.*nyuy(:,1) .* exp(-II.*2.*pi.*nyuy(:,1).*ty(3)).*(ty(3)-ty(1)) ...
               + exp(-II.*2.*pi.*nyuy(:,1).*ty(3)) - exp(-II.*2.*pi.*nyuy(:,1).*ty(1))) ...
              .*(sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)) ...
              ./ (4.*pi.^2.*nyuy(:,1).^2 .* (ty(1)-ty(2)) .* (ty(1)-ty(3))) ;
end

if (abs(ty(1)-ty(3)))<3.*eps
    % this works well, 4/10/2009
    GT(:,nrd) = (tx(1).*(ty(3)-ty(2)) + tx(3).*(ty(2)-ty(1))) ...
              .*(II.*2.*pi.*nyuy(:,1) .* exp(-II.*2.*pi.*nyuy(:,1).*ty(3)).*(ty(2)-ty(3)) ...
               + exp(-II.*2.*pi.*nyuy(:,1).*ty(2)) - exp(-II.*2.*pi.*nyuy(:,1).*ty(3))) ...
              .*(sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)) ...
              ./ (4.*pi.^2.*nyuy(:,1).^2 .* (ty(1)-ty(2)) .* (ty(2)-ty(3))) ;
end

if (abs(ty(1)-ty(2)))<3.*eps
    % this works well, 4/10/2009
    GT(:,nrd) = (tx(1).*(ty(3)-ty(2)) + tx(2).*(ty(1)-ty(3))) ...
              .*(II.*2.*pi.*nyuy(:,1) .* exp(-II.*2.*pi.*nyuy(:,1).*ty(1)).*(ty(1)-ty(3)) ...
               + exp(-II.*2.*pi.*nyuy(:,1).*ty(1)) - exp(-II.*2.*pi.*nyuy(:,1).*ty(3))) ...
              .*(sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)) ...
              ./ (4.*pi.^2.*nyuy(:,1).^2 .* (ty(2)-ty(3)) .* (ty(1)-ty(3))) ;
end

% formula in (f,g) = (0,0)
GT(nrd,nrd) = 1./2 .*(tx(1).*(ty(3)-ty(2)) + tx(2).*(ty(1)-ty(3)) + tx(3).*(ty(2)-ty(1))) ...
            .*(sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)) ;

% unexpected NaN and Inf removal by coordinate rotation

cisin = isnan(GT) ;
cisii = isinf(GT) ;

if sum(cisin(:))+sum(cisii(:)) > 0

    rot = 5.0e-7 ;
    nyuxt = cos(rot).*nyux+sin(rot).*nyuy ;
    nyuyt = -sin(rot).*nyux+cos(rot).*nyuy ;

    G_dom = 4.* pi.^2 ...
          .* (nyuxt .* (tx(1)-tx(2)) + nyuyt .* (ty(1)-ty(2))) ...
          .* (nyuxt .* (tx(3)-tx(1)) + nyuyt .* (ty(3)-ty(1))) ...
          .* (nyuxt .* (tx(3)-tx(2)) + nyuyt .* (ty(3)-ty(2))) ;
    
    G_num = (tx(1).*(ty(3)-ty(2)) + tx(2).*(ty(1)-ty(3)) + tx(3).*(ty(2)-ty(1))) ...
          .*(exp(-II.*2.*pi.*(nyuxt.*tx(1) + nyuyt.*ty(1))) .* (nyuxt.*(tx(3)-tx(2)) + nyuyt.*(ty(3)-ty(2))) ...
           + exp(-II.*2.*pi.*(nyuxt.*tx(2) + nyuyt.*ty(2))) .* (nyuxt.*(tx(1)-tx(3)) + nyuyt.*(ty(1)-ty(3))) ...
           + exp(-II.*2.*pi.*(nyuxt.*tx(3) + nyuyt.*ty(3))) .* (nyuxt.*(tx(2)-tx(1)) + nyuyt.*(ty(2)-ty(1)))) ;

%     GT = (sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)).*G_num./G_dom ;

    GT(cisin) = (sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)).*G_num(cisin)./G_dom(cisin) ;
    GT(cisii) = (sqrt(tra_int).*exp(II.*phase) - sqrt(tra_back).*exp(II.*pha_back)).*G_num(cisii)./G_dom(cisii) ;

end

y = GT ;

% ------------------------------------------------------
% Sinc function
% ------------------------------------------------------
function y = sinc(xx)

y = ones(size(xx)) ;
is = find(xx) ;
y(is) = sin(pi*xx(is))./(pi*xx(is)) ;

% end of this file - Kenji Yamazoe

function polyg = cc_or_ccw(polyg)
% Function M-file cc_or_ccw.m
% *************************************************************************
% Modified for the University of California: 5/9/2010
% Last modified date: 10/2/2014
%
% - Use this code at your own risk.
% - Freely used inside the University of California.
% - Do not use commercially.
% - If you make a significant improvement, please let me know to 
%   kenji.yamazoe@gmail.com.
%
% *************************************************************************
% ****** Copyright (C) Kenji Yamazoe, 2004-2014. All Right Reserved. ******
% *************************************************************************

num_poly = length(polyg) ;
for i_poly = 1:num_poly
    
    [irow,icol] = size(polyg(i_poly).xy) ;
    
    dp_x = [polyg(i_poly).xy(1,1:icol),polyg(i_poly).xy(1,1),polyg(i_poly).xy(1,2)] ;
    dp_y = [polyg(i_poly).xy(2,1:icol),polyg(i_poly).xy(2,1),polyg(i_poly).xy(2,2)] ;
    
    csign = 0 ;
    for i=1:icol
        r1_vec = [dp_x(i+1)-dp_x(i),dp_y(i+1)-dp_y(i),0] ;
        r2_vec = [dp_x(i+2)-dp_x(i+1),dp_y(i+2)-dp_y(i+1),0] ;
        cross_v = cross(r1_vec,r2_vec) ;
        csign = csign + sign(cross_v(3)) ;
    end
    
    if csign>0
        polyg(i_poly).xy = fliplr(polyg(i_poly).xy) ;
    end
    if csign ==0
        fprintf('Please make sure the rotation direction manually: polygon %d.\n',i_poly) ;
    end
    
end

% - end of this file, Kenji Yamazoe



function uieDefocus_Callback(hObject, eventdata, handles)
% hObject    handle to uieDefocus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uieDefocus as text
%        str2double(get(hObject,'String')) returns contents of uieDefocus as a double


% --- Executes during object creation, after setting all properties.
function uieDefocus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uieDefocus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function uieAngle_Callback(hObject, eventdata, handles)
% hObject    handle to uieAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uieAngle as text
%        str2double(get(hObject,'String')) returns contents of uieAngle as a double


% --- Executes during object creation, after setting all properties.
function uieAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uieAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function uieAzimuth_Callback(hObject, eventdata, handles)
% hObject    handle to uieAzimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uieAzimuth as text
%        str2double(get(hObject,'String')) returns contents of uieAzimuth as a double


% --- Executes during object creation, after setting all properties.
function uieAzimuth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uieAzimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in uibPupilMask.
function uibPupilMask_Callback(hObject, eventdata, handles)
% hObject    handle to uibPupilMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mask
pupil = getappdata(gcf,'pupil');
pupilAmp = abs(pupil);
mask = ones(size(pupil));
mask(pupilAmp<max(pupilAmp(:))/10) = 0;
xp_um = getappdata(gcf,'xp_um');
yp_um = getappdata(gcf,'yp_um');
imagesc(handles.pupilAmp,xp_um,yp_um,mask); colorbar(handles.pupilAmp);
xlabel(handles.pupilAmp,'x/um'),ylabel(handles.pupilAmp,'y/um');title(handles.pupilAmp,'Pupil mask');
% Hint: get(hObject,'Value') returns toggle state of uibPupilMask


% --- Executes on button press in uicbRemoveTilt.
function uicbRemoveTilt_Callback(hObject, eventdata, handles)
% hObject    handle to uicbRemoveTilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of uicbRemoveTilt


% --- Executes on button press in uicbRemoveDef.
function uicbRemoveDef_Callback(hObject, eventdata, handles)
% hObject    handle to uicbRemoveDef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of uicbRemoveDef


% --- Executes on button press in uicbUnwrapPhase.
function uicbUnwrapPhase_Callback(hObject, eventdata, handles)
% hObject    handle to uicbUnwrapPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of uicbUnwrapPhase


% --- Executes on button press in uibLoadMask.
function uibLoadMask_Callback(hObject, eventdata, handles)
% hObject    handle to uibLoadMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mask
[fn,pn]=uigetfile({'*.mat','Pupil mask (*.mat)'},'Loading','data\mask');
fileformat=fn(end-2:end);
filename= strcat(pn,fn);
switch fileformat
    case 'mat'
        load(filename);
end
xp_um = getappdata(gcf,'xp_um');
yp_um = getappdata(gcf,'yp_um');
imagesc(handles.pupilAmp,xp_um,yp_um,mask); colorbar(handles.pupilAmp);
xlabel(handles.pupilAmp,'x/um'),ylabel(handles.pupilAmp,'y/um');title(handles.pupilAmp,'Pupil mask');

% --- Executes on button press in uibAnalyze.
function uibAnalyze_Callback(hObject, eventdata, handles)
% hObject    handle to uibAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mask
pupil = getappdata(gcf,'pupil');
pupilAmp = abs(pupil);
pupilPha = getappdata(gcf,'pupilPha');
xp_um = getappdata(gcf,'xp_um');
yp_um = getappdata(gcf,'yp_um');
if length(mask)~=length(pupil)
    mask = ones(length(pupil));
end
if get(handles.uicbUnwrapPhase,'Value')
    pupilPha = GDS.utils.UnwrapPhaseBySortingReliabilityWithMask(pupilPha,mask*255);
end
pupilPha(mask==0) = NaN;
if get(handles.uicbRemoveSph,'Value')
    propdis_um=str2double(get(handles.distance,'String'))*1e3;
    NA = sin(atan((max(xp_um)-min(xp_um))/2/propdis_um));
    pupilPha = GDS.utils.removeSphere(pupilPha,NA,mask);
end
if get(handles.uicbRemoveTilt,'Value')
    pupilPha = GDS.utils.DelTilt(pupilPha);
end
if get(handles.uicbRemoveDef,'Value')
    pupilPha = GDS.utils.DelDefocus(pupilPha);
end

setappdata(gcf,'pupilPhaA',pupilPha);
RMS = std(pupilPha(mask==1))/2/pi;
set(handles.uitRMS,'String',num2str(RMS));
imagesc(handles.pupilAmp,xp_um,yp_um,pupilAmp); colorbar(handles.pupilAmp);
xlabel(handles.pupilAmp,'x/um'),ylabel(handles.pupilAmp,'y/um');title(handles.pupilAmp,'Pupil amplitude');
imagesc(handles.pupilPhase,xp_um,yp_um,pupilPha/2/pi); colorbar(handles.pupilPhase);
xlabel(handles.pupilPhase,'x/um'),ylabel(handles.pupilPhase,'y/um');title(handles.pupilPhase,'Pupil phase');


% --- Executes on button press in uicbRemoveSph.
function uicbRemoveSph_Callback(hObject, eventdata, handles)
% hObject    handle to uicbRemoveSph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of uicbRemoveSph



function uieNZern_Callback(hObject, eventdata, handles)
% hObject    handle to uieNZern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uieNZern as text
%        str2double(get(hObject,'String')) returns contents of uieNZern as a double


% --- Executes during object creation, after setting all properties.
function uieNZern_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uieNZern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in uibDecom.
function uibDecom_Callback(hObject, eventdata, handles)
% hObject    handle to uibDecom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mask
pupilPha = getappdata(gcf,'pupilPhaA');
if isempty(pupilPha)
   pupilPha = getappdata(gcf,'pupilPha'); 
end
NA = str2double(get(handles.uieNA,'String'));
xp_um = getappdata(gcf,'xp_um');
yp_um = getappdata(gcf,'yp_um');
propdis_um=str2double(get(handles.distance,'String'))*1e3;
R_um = propdis_um*tan(asin(NA));
[x,y] = meshgrid(xp_um,yp_um);
[sx,sy] = size(pupilPha);
order =str2double(get(handles.uieNZern,'String'));
X=x/R_um;
Y=y/R_um;
basis = zeros(sx*sy, order + 1);
basis(:,1) = mask(:);
for k = 1:order
    [n, m] = j2nm(k);
    zRTH = ZgenNM(n,m);
    [TH, R] = cart2pol(X(:), Y(:));
    basis(:, k+1) = zRTH(R,TH).*mask(:);
end
dZrn =pinv(basis(mask==1,:))*pupilPha(mask==1);
figure(100),bar(0:order,dZrn);
xlabel('Zernike terms'),ylabel('Coefficients / Waves');
% set(h,'LineWidth',2); 
set(gca,'FontSize',14);



function uieNA_Callback(hObject, eventdata, handles)
% hObject    handle to uieNA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uieNA as text
%        str2double(get(hObject,'String')) returns contents of uieNA as a double


% --- Executes during object creation, after setting all properties.
function uieNA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uieNA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Z = ZgenNM(n,m)
R = RgenNM(n,abs(m));
A = AgenNM(m);
Z = @(r, phi) R(r).*A(phi);

function [n,m] = j2nm(j)
% j = j+1; % 0 is piston
smct = 2;
ct = 0;
numInSmct = 3;
while j > 0
    smct = smct + 2;
    j = j - numInSmct;
    numInSmct = smct + 1;
    ct = ct + 1;
end
n = 2*ct;
m = 0;

for k = 1:abs(j)
    if isodd(k)
        n = n - 1;
        m = m + 1;
    end
end
if isodd(abs(j))
    m = -m;
end

function f = RgenNM(n,m)
p = (n-m)/2;
f = @(r) 0;

for k = 0:p
    coef = (-1)^k*nchoosek(n-k,k)*nchoosek(n-2*k,(n-m)/2-k);
    ex = (n - 2*k);
    f = @(r) f(r) + coef*r.^ex;
end


% Forms azimuthal function based on n,m
function g = AgenNM(m)
if m > 0
    g = @(phi) cos(m*phi);
elseif m < 0
    g = @(phi) -sin(m*phi);
else
    g = @(phi) 1;
end


% Generates a coefficient vector
function coefs = RgenNMCoef(n,m)
p = (n-m)/2;

coefs = zeros(1,p + 1);
for k = 0:p
    coef = (-1)^k*nchoosek(n-k,k)*nchoosek(n-2*k,(n-m)/2-k);
    ex = (n - 2*k);
    coefs(ex + 1) = coef;
end

% Generates a coefficient vector
function angCoefs = AgenNMCoef(m)
cosCoefs = zeros(1,abs(m));
sinCoefs = zeros(1,abs(m));
if m > 0
    cosCoefs(m) = 1;
elseif m < 0
    sinCoefs(abs(m)) = -1;

elseif m == 0
    angCoefs = 1;
    return
end
angCoefs = [cosCoefs; sinCoefs];
angCoefs = [0;angCoefs(:)];


% --- Executes on selection change in uipPropMethod.
function uipPropMethod_Callback(hObject, eventdata, handles)
% hObject    handle to uipPropMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns uipPropMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from uipPropMethod


% --- Executes during object creation, after setting all properties.
function uipPropMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipPropMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
