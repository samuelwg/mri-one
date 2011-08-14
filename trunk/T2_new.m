function [ M0, T2_res, R_square percentage ] = T2_new( file_name, CPMG )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T2 = load (file_name);
sol_perc_num = filename_parser (file_name);

if (CPMG == 1)
    exper = 'CPMG';
    jump  = 2;
else
    exper = 'CP';
    jump  = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%T2 Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%
% we'll find the zero-crossing of differential from positive to negative

[C I] = max (T2(:,2));
max_Voltage = C;
t_delay = T2(I, 1);
figure;
plot(T2(:,1), T2(:, 2));
deriv_Volt = diff(T2(:, 2))./diff(T2(:,1)); % we obtained the differential
deriv_Volt = [0; deriv_Volt];
%% should add afterwards
%% we'll detect zero-crossing:
for i = 
end

