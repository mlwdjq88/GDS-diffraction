function p_amp = sphereProp(polyg,x_nm,y_nm,z_nm,lambda_nm,offsetAngle,azimuth)
p_amp =zeros(size(x_nm));
num_poly = length(polyg);
parfor_progress(num_poly); % progress bar
try
    parfor i_poly = 1:num_poly
        parfor_progress; % get progress
%         if mod(i_poly,10)==1
%             fprintf('%d/%d propagation is done\n',i_poly,num_poly);
%         end
        xy0 = mean(polyg(i_poly).xy,2);
        x = xy0(1);
        y = xy0(2);
        z =x*tan(offsetAngle)*cos(azimuth)+y*tan(offsetAngle)*sin(azimuth);
        R_nm = sqrt((x-x_nm).^2+(y-y_nm).^2 +(z-z_nm).^2);
        p_amp = p_amp + polyg(i_poly).tx.*exp(1i*polyg(i_poly).phase+1i*2*pi/lambda_nm*R_nm)./R_nm.^2*(z_nm-z);
    end
catch
    for i_poly = 1:num_poly
        parfor_progress; % get progress
%         if mod(i_poly,10)==1
%             fprintf('%d/%d propagation is done\n',i_poly,num_poly);
%         end
        xy0 = mean(polyg(i_poly).xy,2);
        x = xy0(1);
        y = xy0(2);
        z =x*tan(offsetAngle)*cos(azimuth)+y*tan(offsetAngle)*sin(azimuth);
        R_nm = sqrt((x-x_nm).^2+(y-y_nm).^2 +(z-z_nm).^2);
        p_amp = p_amp + polyg(i_poly).tx.*exp(1i*polyg(i_poly).phase+1i*2*pi/lambda_nm*R_nm)./R_nm.^2*(z_nm-z);
    end
end
parfor_progress(0);