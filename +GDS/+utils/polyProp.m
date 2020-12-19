function p_amp = polyProp(polyg,backg,nrd,nyux,nyuy)
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
p_amp =zeros(size(nyux));

% The number of the polygons
num_poly = length(polyg) ;
xunit_vec = [1,0] ;
parfor_progress(num_poly); % progress bar

try
    parfor i_poly = 1:num_poly
        parfor_progress; % get progress
%         if mod(i_poly,10)==1
%             fprintf('%d/%d propagation is done\n',i_poly,num_poly);
%         end
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
        end
    end
catch
    for i_poly = 1:num_poly
        parfor_progress; % get progress
%         if mod(i_poly,10)==1
%             fprintf('%d/%d propagation is done\n',i_poly,num_poly);
%         end
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
        end
    end
end
parfor_progress(0);
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