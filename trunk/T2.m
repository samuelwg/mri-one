function [ M0, T2_res, R_square percentage ] = T2( file_name, CPMG)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
const_diff = 20;

T2 = load (file_name);
sol_perc_num = filename_parser (file_name);

if (CPMG == 1)
    exper = 'CPMG';
    jump  = 2;
else
    exper = 'CP';
    jump  = 1;
end

%%%% looking for t=0
t0_point = find(T2(:, 1) >0, 1 );
t0_point_ind = t0_point(1)
T2 = T2(t0_point_ind:end, :);

%%% if this is 1% solution, the final point is 0.06
if (sol_perc_num == 1)
    tfin_point = find(T2(:, 1) <0.06, 1 );
    tfin_point_ind = t0_point(end);
    T2 = T2(1:tfin_point_ind, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% T_2 Calculation           %%%%%%%%%%%%%%

 %%% Let us take a look on T1 data (CP)
[C I] = max (T2(:,2));
max_Voltage = C;
t_delay = T2(I, 1);
figure;
plot(T2(:,1), T2(:, 2));
title(file_name);
hold on;
deriv_Volt = diff(T2(:, 2))./diff(T2(:,1)); % we obtained the differential
deriv_Volt = [0; deriv_Volt];
Threshold = find_thresh(deriv_Volt);


% for the envelope we will look for the maximum between two points that 
% the derivative is above the threshold
deriv_threshed = (deriv_Volt>Threshold);
d_deriv_border = find(diff (deriv_threshed)>0)
i1 = 1;
for i = 1:length(d_deriv_border)-1
    if d_deriv_border(i+1)-d_deriv_border(i)>const_diff
        d_deriv_border_new(i1) = d_deriv_border(i);
        i1= i1 +1;
    end
end
d_deriv_border = d_deriv_border_new;
plot(T2(d_deriv_border,1),T2(d_deriv_border, 2), 'ro');

figure;
plot(T2(:,1), deriv_Volt);
j = 1;
for i=1:jump:length(d_deriv_border)-1
    [C I] = max (T2(d_deriv_border(i):d_deriv_border(i+1), 2));
    T_2_res_y(j) = C;
    T_2_res_x(j) = T2(I+d_deriv_border(i), 1);
    j = j + 1;
end
figure;
plot (T_2_res_x, T_2_res_y);
title (['T_2 calculation out of ',num2str(sol_perc_num),'% ', exper,' experiment']);
legend ([exper, ' experimental data']);
% hold on
% 
% s2 = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0, 0],'Upper', [Inf, Inf], 'Startpoint',[1e-6, max_Voltage]);
% h2 = fittype ('M0*exp(-(t-t_del)/T_2)', 'coefficients', {'T_2', 'M0'}, 'independent', 't', 'problem', 't_del', 'options', s2)
% [c2, gof2] = fit (T_2_res_x', T_2_res_y', h2, 'problem', t_delay)
% 
% plot(c2, 'k');
% xlabel ('t [sec]');
% ylabel ('M_{trans} [Volt]');
% text ('Position', [3e-3, 0.4], 'String', {strcat(' Fit results: M_0=', num2str(c2.M0)), strcat(' {T_2}=', num2str(c2.T_2))});

% for a comparison we'll try to perform fit for the linear function after
% logarithm evauation
figure;
plot (T_2_res_x, log(T_2_res_y), 'o');
title (['T_2 calculation out of ', num2str(sol_perc_num),'% ', exper,'experiment (linear fitting curve)'])
legend ('experimental data')
hold on

%s4 = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0, 0],'Upper', [Inf, Inf], 'Startpoint',[1e-6, max_Voltage]);
h5 = fittype ('poly1');
[c5 gof5] = fit (T_2_res_x', log(T_2_res_y)', h5)
plot(c5, 'k');
xlabel ('t [sec]');
ylabel ('ln(M_{trans}) [Volt]');
text ('String', {strcat(' Fit results: M_0=', num2str(exp(c5.p2))), strcat(' {T_2}=', num2str(-1/c5.p1))});

M0 = exp(c5.p2);
T2_res = -1/c5.p1;
R_square = gof5.rsquare;
percentage = sol_perc_num;


end

