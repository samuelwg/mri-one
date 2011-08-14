function [ M0, T2_res, R_square percentage ] = T2_new( file_name, CPMG )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
neib_points = 3;
T2 = load (file_name);
sol_perc_num = filename_parser (file_name);

if (CPMG == 1)
    exper = 'CPMG';
    jump  = 2;
else
    exper = 'CP';
    jump  = 1;
end

%%%% we will find here t=0 point
t0_point = find(T2(:, 1) >0, 1 );
t0_point_ind = t0_point(1)


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
zerocross_list = [];
j=1
for i = t0_point_ind:length(T2(:,2))
    if(sign(deriv_Volt(i)) > sign(deriv_Volt(i+1)))
        zerocross_list(j) = i+t0_point_ind %we are saving the 'before' point
    end
end
%now we will look for the maximum in several (neib_points) points
local_maxima_x = [];
local_maxima_y = [];
l = 1;
for k = 1:jump:length(zerocross_list)
    [C I] = max(T2(zerocross_list(k)-neib_points:zeroicross_list(k)+neib_points, 2));
    local_maxima_y (l) = C;
    local_maxima_x (l) = T2(I+zerocross_list(k), 1);
    l = l+1;
end
figure;
plot(local_maxima_x, local_maxima_y);
