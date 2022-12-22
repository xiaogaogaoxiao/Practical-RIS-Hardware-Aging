function [S_mn]=function_Smn(p_BS, centers_IRS, p_vehicle_trajectory, total_move_time, dx, dy, isotropic, exp_phi) 
% Author: Ke(Ken)WANG from Macao Polytechnic University
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Update infomation: v0.1(2020/11/20), v0.4(2022/04/18)
%
% This function aims to calculate S_mn(t), i.e., eq(8) in the GC paper 2021
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% original article.
%
% Example:
%
% TBD

    global c0 fc;
    
    lambda_c = c0/fc;
    
    %This gain is $G_T^{mn}$
    %This is the gain from BS to IRS, that means this gain always euqals to 1
    %since we assume the transceiver is isotropic
    G_BS_to_IRS_direction = function_antenna_gain_TR(p_BS, centers_IRS);
    
    %This gain is $G_{mn}^T$. If the IRS is isotropic, the gain is 1, otherwise we use
    %"function_gain_IRS" to calculate the gain. 
    G_IRS_to_BS_direction = function_gain_IRS(centers_IRS , p_BS, dx, dy, isotropic);
    
    %Total gain between BS and IRS
    A_BS_IRS = lambda_c/(4*pi) * sqrt(G_BS_to_IRS_direction.*G_IRS_to_BS_direction)./vecnorm((centers_IRS-p_BS).'); 

    S_mn = zeros(1, total_move_time);

    G_Vehicle_to_IRS_direction=[];

    G_IRS_to_Vehicle_direction=[];

    %This for-loop aims to calculate $G_R^{mn}(t)$ and $G_{mn}^R(t)$ in eq(6)
    %and its above content in the GC paper 2021
    for tt = 1 : total_move_time

        %This is the gain from Vehicle to IRS, $G_R^{mn}(t)$. That means this gain always euqals to 1
        %since we assume the transceiver is isotropic
        temp_G_Vehicle_to_IRS_direction = function_antenna_gain_TR(p_vehicle_trajectory(tt,:), centers_IRS); 
        G_Vehicle_to_IRS_direction =[G_Vehicle_to_IRS_direction; temp_G_Vehicle_to_IRS_direction];

        %This is the gain from IRS to Vehicle, $G_{mn}^R(t)$. If the IRS is isotropic, the gain is 1, otherwise we use
        %"function_gain_IRS" to calculate the gain. 
        temp_G_IRS_to_Vehicle_direction = function_gain_IRS(centers_IRS , p_vehicle_trajectory(tt,:), dx, dy, isotropic); 
        G_IRS_to_Vehicle_direction =[G_IRS_to_Vehicle_direction; temp_G_IRS_to_Vehicle_direction];

    end
    
    %This for-loop aims to calculate S_mn(t), i.e., eq(8) in the GC paper 2021
    %Note that the $\sqrt{P_t}$ is not included here. Instead, we multiply
    %it when we finally calculate the received power 
    for tt = 1 : total_move_time
    
        a_IRS_Vehicle = ...
            lambda_c/(4*pi)*sqrt(G_Vehicle_to_IRS_direction(tt,:).*G_IRS_to_Vehicle_direction(tt,:))./ ...
            vecnorm((p_vehicle_trajectory(tt,:) - centers_IRS).');
        
        %The phase shift by IRS, i.e., "exp_phi(:,t)", is given by main.m
        %For example, "exp_phi_optimal = exp(-1i*2*pi*fc.*tau_mn-1i.*phi_optimal)"
        temp_S_mn = sum(A_BS_IRS .* a_IRS_Vehicle .* exp_phi(:,tt).');
        
        %eq(8) in the GC paper 2021
        S_mn(tt) = temp_S_mn;
        
        clear temp_S_mn

    end

end
