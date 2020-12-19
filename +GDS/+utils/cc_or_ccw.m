
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
