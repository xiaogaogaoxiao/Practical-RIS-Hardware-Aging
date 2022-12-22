function  [exp_phi_out] = function_hardware_aging(beta_exp_phi_wo_HA, total_move_time, N, failure_rate_each_time, a, b, c)
% Author: Ke(Ken)WANG from Macao Polytechnic University
% Email: ke.wang@mpu.edu.mo, kewang0225@gmail.com
% Update infomation: v0.2(2022/Dec)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% original article.
%
% Example: TBD
HA_vector_1_each_time = randsrc(N, 1, [0, 1; failure_rate_each_time, 1-failure_rate_each_time]);

% size(HA_vector_1) = N x total_move_time 
% HA_vector_1 = repmat(HA_vector_1_each_time, 1, total_move_time);

HA_0_position_each_time = find(HA_vector_1_each_time==0); 

HA_vector_only_01 = repmat(HA_vector_1_each_time, 1, total_move_time);

for ii = 1 : length(HA_0_position_each_time)
    temp_0 = random('Uniform',0,2*pi);
    temp_1 = random('Uniform',b,1);
    HA_vector_1_each_time(HA_0_position_each_time(ii)) = temp_1.*exp(-1i*temp_0); % 幅度范围现在为从b到1
    clear temp_0 temp_1
end

HA_vector_1 = repmat(HA_vector_1_each_time, 1, total_move_time);
HA_vector_2 = HA_vector_1;
HA_vector_3 = HA_vector_1;
HA_vector_3(HA_vector_3==1)=0;

exp_phi_out = [];

exp_phi_out = HA_vector_only_01 .* beta_exp_phi_wo_HA + HA_vector_3; 

end