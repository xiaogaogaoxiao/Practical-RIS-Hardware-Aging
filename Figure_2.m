% Author: Ke(Ken)WANG from Macao Polytechnic University
% Email: ke.wang@mpu.edu.mo, kewang0225@gmail.com
% Update infomation: v1.0(Dec/2022)


%% Clean All & Timer Begins
close all;
clear;
tic;

%% System Parameters Initialization

global c0 fc speed; % Global parameters

% You CAN change the parameters below to obtain different figures
% In this simulation we IGNORE the phase drifts

kappa_t = 0; % Transmitter HWI proportionality coefficient 

kappa_r = kappa_t; % Receiver HWI proportionality coefficient

q = 2;
alpha_1 = -pi*(2^(-q)); % RIS HWI, gamma ~ U (alpha_1, alpha_2) % save life_time_400_5000 'SE_all_different_failure_rate' save SE_direct 'SE_direct'
alpha_2 = pi*(2^(-q));
alpha = (abs(alpha_1)+abs(alpha_2))/2;

realization_HWI = 1; % Number of realizations of HWI. In fact, realization_HWI = 1 is enough since the element number >> 0 

sqrt_N =  64; % IRS has sqrt_N x sqrt_N elements

bias_x_axis = 0; % Bias of IRS on the x-axis, normally it is set to 0

c0 = physconst('Lightspeed'); % Light speed

fc= 2.4e9; % Carrier frequency is 2.4 GHz

disp(['Step 1: Parameter Initialization begin. In this setup, the IRS has ' num2str(sqrt_N) ' x ' num2str(sqrt_N) ' elements.']);

dx = (c0/fc)/(2*sqrt(pi)); dy = dx;

h_IRS = 10; % Height of the IRS

h_User = 1.5; % Height of the vehicle

h_BS = 15; % Height of the BS

speed = 0; % Speed that the user moves

total_time_slot = 1; % Total time for the moving by the user

Pt = 0.1; % Transmit power is 0.1W = 20dBm

sigma_square = 10e-12; % AWGN noise, sigma^2 = -80dBm

% The coordinate of the user
x_User_start = 0;
y_User = 0+h_User;
z_User = 50; 

% The coordinate of the Base Station, i.e., the transmitter
x_BS = -50;
y_BS = 0+h_BS;
z_BS = 20; 

% -------- Hardware aging parameters begin -------- 

a=1; % In this paper, a == 1;

b=1;  % 0~1, practical value is 0.2

c=0.43*pi;  % practical value is 0.43pi

total_running_time_t = 5*1e3; % hours

life_time_L = 2*1e2; % hours 

max_shape_para_rho = 3.5;

early_external_failures = 0; 

realization_HWI_HA = 5e3; % Number of realizations of HA, you can set 1e3

% -------- Hardware aging parameters end -------- 

%% Compute the trajectory for one vehicle pass

disp(['Step 2: Compute the trajectory for one vehicle pass. The speed is ' num2str(speed) ' m/s.']);

p_BS = [x_BS, y_BS, z_BS];

p_User_start = [x_User_start, y_User, z_User];

% Total trajectory = starting point + moving trajectory per second (we start from 1s rather than 0s)
p_User_trajectory = [];

% The vehicle starts at p_vehicle_start, and travels at the speed during
% total_time period
for t_moves = 1 : total_time_slot
    
    temp_p_User_trajectory =  p_User_start+function_User_moving_xdir(speed, t_moves);
    
    p_User_trajectory = [p_User_trajectory; temp_p_User_trajectory];
    
end

% 'A_0' is the gain for the direct path. It is a theoretical result.
[A_0] = function_A0(p_BS, p_User_trajectory);

%% Compute some basic results

disp('Step 3: Compute some basic results.');

% 1/0 stands for isotropic IRS/practical IRS. In this paper, we only
% consider isotropic case. 
isotropic = 1;

% Locate the element positions for IRS
% size(centers_IRS) = sqrt_N*sqrt_N x 3
[centers_IRS] = function_centers_IRS(sqrt_N, sqrt_N, dx, dy);

% Add the height of IRS
centers_IRS(:, 2) = centers_IRS(:, 2) + h_IRS;

% Add the bias of IRS on the x-axis
centers_IRS(:, 1) = centers_IRS(:, 1) + bias_x_axis;

% 'A_mn_theoretical' is the gain for the IRS path, which is a theoretical result.
[A_mn, A_mn_matrix] = function_Amn(p_BS, centers_IRS, p_User_trajectory, total_time_slot, dx, dy, isotropic);

% Compute the delays, size(tau_0) = 1 x 50
[tau_0, tau_n] = function_time_delay(centers_IRS, p_BS, p_User_trajectory); 

% Optimal phase shift
phi_optimal = 2*pi*(fc*(tau_0-tau_n)+ceil(-(fc*(tau_0-tau_n))));

[signal_direct_HWI]=function_S0(A_0, tau_0);

%% Compute the gain when the IRS is non-isotropic, after phase optimization, w/ HWI and w/ HA

disp('Step 4: Compute the gain when the IRS is non-isotropic, after phase optimization, w/ HWI and w/ HA.');

% size(failure_rate) = 1 x total_running_time_t
failure_rate =function_failure_rate(total_running_time_t, life_time_L, max_shape_para_rho, early_external_failures); 

failure_rate(failure_rate>1) = 1;

index_failure_rate_vector = linspace(1, total_running_time_t, 10);

for kk = 1 : length(ceil(index_failure_rate_vector))
    aa = ceil(index_failure_rate_vector);
    bb = aa(kk);
    failure_rate_vector(kk) = failure_rate(bb);
end

SE_all_different_failure_rate = zeros(1, length(failure_rate_vector));

for ii = 1 : length(failure_rate_vector)

% for-loop for realizations of HA
for jj = 1 : realization_HWI_HA

disp(['           =============This is the ',num2str(jj), '-th realiztion of HWI&HA=============']);

% Generate Theta_IRS_HWI
% RIS HWI of each element, size(gamma_n) = sqrt_N*sqrt_N x total_move_time
gamma_n = random('Uniform', alpha_1, alpha_2, sqrt_N*sqrt_N, total_time_slot);

%disp(['           =====The Theta_IRS_HWI w/ HA of the (1, 1)-th element is ',num2str(gamma_n(1,1)), '=====']);

phi_IRS_HWI = -1i*2*pi*fc*tau_n - 1i*(phi_optimal+gamma_n);  

exp_phi_IRS_HWI = exp(phi_IRS_HWI);

beta_n_IRS_HWI = function_beta_n((phi_optimal+gamma_n), a, b, c);

beta_exp_phi_IRS_HWI = beta_n_IRS_HWI.*exp_phi_IRS_HWI;

%size(beta_exp_phi_IRS_HWI_HA) = N x total_move_time
[beta_exp_phi_IRS_HWI_HA] = function_hardware_aging(beta_exp_phi_IRS_HWI,total_time_slot, sqrt_N*sqrt_N, failure_rate_vector(ii),a,b,c);

% 'signal_IRS_HWI' is the gain with HWI for cascaded link, it is a Monte-Carlo simulation result
[signal_IRS_HWI_HA] = function_Smn(p_BS, centers_IRS, p_User_trajectory, total_time_slot, dx, dy, isotropic, beta_exp_phi_IRS_HWI_HA);

% Total gain with HWI, it is a Monte-Carlo result
% size(signal_HWI_HA_total) = 1 x total_move_time
signal_HWI_HA_total = signal_IRS_HWI_HA + signal_direct_HWI;

signal_HWI_HA_total_all_Realization(jj, :) = signal_HWI_HA_total;

clear signal_IRS_HWI_HA signal_HWI_HA_total

end

[SE_IRS_HWI_HA_MC] = function_compute_SE(signal_HWI_HA_total_all_Realization, Pt, sigma_square, kappa_t, kappa_r);

SE_all_different_failure_rate(ii) = SE_IRS_HWI_HA_MC;

clear SE_IRS_HWI_HA_MC

end

[SE_DIRECT] = function_compute_SE(signal_direct_HWI, Pt, sigma_square, kappa_t, kappa_r);

SE_direct = ones(1,size(index_failure_rate_vector,2)) * SE_DIRECT;

%% Plot the simulation results

close all;

disp('Step 8: Plot the simulation results.');

% Plot
figure(1); hold on; box on; grid on;

total_running_time_t_vector = index_failure_rate_vector;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 

plot(total_running_time_t_vector, SE_all_different_failure_rate,'r-^', 'LineWidth', 1.5, 'MarkerSize', 10);
plot(total_running_time_t_vector, SE_direct, 'k-s', 'LineWidth', 1.5, 'MarkerSize', 10);

legend('$L  = 200$ Hours', 'w/o RIS', 'Interpreter','LaTex','Location','NorthEast');

xlabel('Total Runtime $t$ (Hours)','Interpreter','LaTex');
ylabel('Spectral Efficiency (bit/s/Hz)','Interpreter','LaTex');

%% Timer Ends
toc;





