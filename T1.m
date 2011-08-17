function [ M0, T1_res, R_square, percentage] = T1( file_name, iterator )
%   The function calculates 
%   Detailed explanation goes here

sol_perc_num = filename_parser (file_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



T1 = load(file_name);
subplot(2,3, iterator);
plot(T1(:,2),T1(:,1), 'r*');
title (['T_1 experiment for ', num2str(sol_perc_num),'% solution: Inversion recovery']);
legend ('Experimental data')
hold on
%% We'll find an approximate M0 and T1 values:

[C I] = max(T1(:,1));
M0_app = C;
[C I] = min(T1(:,1));
T1_app = T1(I,2)/log(2);

%% Now we'll fit data for the model: 
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf, Inf],...
               'Startpoint',[M0_app T1_app]);
h = fittype('abs(M0*(1-2*exp(-x/T1)))', 'coefficients', {'M0', 'T1'}, 'independent', 'x', 'options', s)
[c, gof] = fit (T1(:,2),T1(:,1), h)
plot(c, 'k');
xlabel ('T_I [msec]');
ylabel ('M_{trans}(T_I) [msec]');
text (...%'Position', [10, 200], 
    'String', strcat(' Fit results: M_0=', num2str(c.M0), ' T_1=', num2str(c.T1)))
M0 = c.M0;
T1_res = c.T1;
R_square = gof.rsquare;
percentage = sol_perc_num;
end

