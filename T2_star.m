function [ M0, T2_res, R_square percentage ] = T2_star( file_name, iterator )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
T2_star = load (file_name);

sol_perc_num = filename_parser (file_name);


[C I] = max (T2_star(:,2));
max_Voltage = C;
ind_start = I + 50; % it is the moment we are starting to fit
t_delay = T2_star(ind_start, 1);
subplot(2, 3, iterator);
plot(T2_star(ind_start:ind_start+820,1), ((T2_star(ind_start:ind_start+820, 2)- T2_star(end, 2))));
legend('experimental data');
hold on

signal_fft = fft(T2_star(:,2))/(length(T2_star(:,2))-ind_start);

s1 = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0, 0],'Upper', [Inf, Inf], 'Startpoint',[1e-6, max_Voltage-T2_star(end, 2)]);

h1 = fittype ('M0*exp(-(t-t_del)/T_2_star)', 'coefficients', {'T_2_star', 'M0'}, 'independent', 't', 'problem', 't_del', 'options', s1)
%[c1, gof1] = fit (T2_star_8percent(ind_start+20:end,1), T2_star_8percent(ind_start+20:end,2), h1)
[c1, gof1] = fit (T2_star(ind_start:ind_start+820,1), T2_star(ind_start:ind_start+820,2)-T2_star(end, 2), h1, 'problem', t_delay)
plot (c1, 'k');
xlabel ('t [sec]');
ylabel ('M_{trans} [Volt]');
%text ('String', {strcat(' Fit results: M_0=', num2str(c1.M0)), strcat(' {T_2}^*=', num2str(c1.T_2_star))});
title (['{T_2}^* experiment for ', num2str(sol_perc_num),'% solution: FID']);

M0 = c1.M0;
T2_res = c1.T_2_star;
R_square = gof1.rsquare;
percentage = sol_perc_num;

end

