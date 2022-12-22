function [M_new, N_new, dx_new, dy_new]=function_parameters_IRS(M, N, dx, dy, isotropic)
% Author: Ke(Ken)WANG from Macao Polytechnic Institute
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Update infomation: v0.1(2020/11/19), v0.2(2021/08/26), v0.3(2021/09/04)
%
% This function aims to reconfigure M, N, dx and dy of an IRS if we set it is
% isotropic. Otherwise those 4 parametres are the same as the INPUT.
% Note that isotropic=1 stands for isotropic IRS, isotropic=0 stands for normal IRS.
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% original article.
%
% Example:
%
% [M_new, N_new, dx_new, dy_new]=function_parameters_IRS(200, 100, 0.1, 0.6)
%
% Please fill in the input parameters.
%
% [M_new, N_new, dx_new, dy_new]=function_parameters_IRS(200, 100, 0.1, 0.6, 3)
%
% Please select the type of IRS. Note that 1 stands for isotropic IRS and 0 stands for normal IRS.
%
% [M_new, N_new, dx_new, dy_new]=function_parameters_IRS(200, 100, 0.1, 0.6, 1)
%
% M_new = 1159
%
% N_new = 579
%
% dx_new = 0.0173
%
% dy_new = 0.1036
%
% [M_new, N_new, dx_new, dy_new]=function_parameters_IRS(200, 100, 0.1, 0.6, 0)
%
% M_new = 200
%
% N_new = 100
%
% dx_new = 0.1
%
% dy_new = 0.6


global c0 fc

lambda_c = c0/fc;

if nargin<5
    
    error('Please fill in the input parameters.')
    
else
    
    %Calculate the area of practical IRS
    area_IRS = M*N*dx*dy;
    
    %Calculate the area of isotropic IRS, note that it is different from the area of practical IRS
    area_isotropic_IRS = lambda_c^2/(4*pi);
    
    %Aspect
    aspect = M/N;
    
    %Recalculate the size of isotropic IRS
    if isotropic == 1
        
        temp_N_new = sqrt(area_IRS/area_isotropic_IRS/aspect);
        M_new = fix(roundn(aspect*temp_N_new, 0));
        N_new = fix(roundn(temp_N_new, 0));
        
        dx_new = (M*dx)/M_new;
        dy_new = (N*dy)/N_new;
    
    %If we don't need isotropic IRS, then we just make the output and the
    %input the same
    elseif isotropic == 0
        
        M_new = M;
        N_new = N;
        dx_new = dx;
        dy_new = dy;
        
    else
        
        error('Please select the type of IRS. Note that 1 stands for isotropic IRS and 0 stands for normal IRS. ')
        
    end
    
end

end
