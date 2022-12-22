function beta_n = function_beta_n(theta_n, a, b, c)
% Author: Ke(Ken)WANG from Macao Polytechnic University
% Email: ke.wang@mpu.edu.mo, kewang0225@gmail.com
% Update infomation: v0.2(2022/Dec)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% original article.
%
% Example: TBD

beta_n = (1-b).*((sin(theta_n-c)+1)/2).^a + b;

end