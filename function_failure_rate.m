function [failure_rate]=function_failure_rate(total_running_time_t, life_time_L, max_shape_para_rho, external_failures)

%temp_shape_para_rho_vector = [ones(1,total_running_time_t*0.1)*1,linspace(1,max_shape_para_rho,total_running_time_t*0.9)];
temp_shape_para_rho_vector = linspace(1,max_shape_para_rho,total_running_time_t);

failure_rate = (temp_shape_para_rho_vector./(life_time_L.^temp_shape_para_rho_vector))...
    .* total_running_time_t.^(temp_shape_para_rho_vector-1)+external_failures;

end

% example 
% a = function_failure_rate(4000, 1000, 3.5, 0.05); plot(a);