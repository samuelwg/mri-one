%%%% NMR Experiment %%%
clear all;
close all;
clc;
%%% loading data for 8 percent solution %%%
Threshold = 9000;
T1_8percent = load('./T1_8percent.csv');
T2_8percent_CPMG = load('./TT2_8percent_CPMG.csv');
T2_8percent_CP = load('./TT2_8percent_CP.csv');
T2_star_8percent = load('./T2_star_8percent.csv');
figure;

%let us calculate the sampling frequency
delta_t = T2_star_8percent(2, 1)-T2_star_8percent(1, 1);
sampling_freq = 1/delta_t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% T_1 Calculation %%%%%%%%%%%%%%%%%%%%%%%


plot(T1_8percent(:,2),T1_8percent(:,1), 'r*');
title ('T_1 experiment for 8% solution: Inversion recovery');
legend ('Experimental data')
hold on
%% We'll find an approximate M0 and T1 values:

[C I] = max(T1_8percent(:,1));
M0_app = C;
[C I] = min(T1_8percent(:,1));
T1_app = T1_8percent(I,2)/log(2);

%% Now we'll fit data for the model: 
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf, Inf],...
               'Startpoint',[M0_app T1_app]);
h = fittype('abs(M0*(1-2*exp(-x/T1)))', 'coefficients', {'M0', 'T1'}, 'independent', 'x', 'options', s)
[c, gof] = fit (T1_8percent(:,2),T1_8percent(:,1), h)
plot(c, 'k');
xlabel ('T_I [msec]');
ylabel ('M_{trans}(T_I) [msec]');
text ('Position', [10, 200], 'String', strcat(' Fit results: M_0=', num2str(c.M0), ' T_1=', num2str(c.T1)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% T_2_star calculation %%%%%%%%%%%%%%%%%

%% We will take a look on the FID:

% figure;
% plot(T2_star_8percent(:,1), T2_star_8percent(:, 2));
%% We will choose the data shortly after the peak
[C I] = max (T2_star_8percent(:,2));
max_Voltage = C;
ind_start = I + 50; % it is the moment we are starting to fit
t_delay = T2_star_8percent(ind_start, 1);
figure;
plot(T2_star_8percent(ind_start:ind_start+820,1), ((T2_star_8percent(ind_start:ind_start+820, 2)- T2_star_8percent(end, 2))));
legend('experimental data');
hold on

signal_fft = fft(T2_star_8percent(:,2))/(length(T2_star_8percent(:,2))-ind_start);

s1 = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0, 0],'Upper', [Inf, Inf], 'Startpoint',[1e-6, max_Voltage-T2_star_8percent(end, 2)]);

h1 = fittype ('M0*exp(-(t-t_del)/T_2_star)', 'coefficients', {'T_2_star', 'M0'}, 'independent', 't', 'problem', 't_del', 'options', s1)
%[c1, gof1] = fit (T2_star_8percent(ind_start+20:end,1), T2_star_8percent(ind_start+20:end,2), h1)
[c1, gof1] = fit (T2_star_8percent(ind_start:ind_start+820,1), T2_star_8percent(ind_start:ind_start+820,2)-T2_star_8percent(end, 2), h1, 'problem', t_delay)
plot (c1, 'k');
xlabel ('t [sec]');
ylabel ('M_{trans} [Volt]');
text ('Position', [4e-4, 0.4], 'String', {strcat(' Fit results: M_0=', num2str(c1.M0)), strcat(' {T_2}^*=', num2str(c1.T_2_star))});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% T_2 Calculation out of CP %%%%%%%%%%%%%%

 %%% Let us take a look on T1 data (CP)
[C I] = max (T2_8percent_CP(:,2));
max_Voltage = C;
t_delay = T2_8percent_CP(I, 1);
figure;
plot(T2_8percent_CP(:,1), T2_8percent_CP(:, 2));
deriv_Volt = diff(T2_8percent_CP(:, 2))./diff(T2_8percent_CP(:,1)); % we obtained the differential
deriv_Volt = [0; deriv_Volt];
figure;
plot(T2_8percent_CP(:,1), deriv_Volt);
% for the envelope we will look for the maximum between two points that 
% the derivative is above the threshold
deriv_threshed = (deriv_Volt>Threshold);
d_deriv_border = find(diff (deriv_threshed)>0);

for i=1:length(d_deriv_border)-1
    [C I] = max (T2_8percent_CP(d_deriv_border(i):d_deriv_border(i+1), 2));
    T_2_results_y(i) = C;
    T_2_results_x(i) = T2_8percent_CP(I+d_deriv_border(i), 1);
end
figure;
plot (T_2_results_x, T_2_results_y);
title ('T_2 calculation out of 8% CP experiment');
legend ('CP experimental data');
hold on

s2 = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0, 0],'Upper', [Inf, Inf], 'Startpoint',[1e-6, max_Voltage]);
h2 = fittype ('M0*exp(-(t-t_del)/T_2)', 'coefficients', {'T_2', 'M0'}, 'independent', 't', 'problem', 't_del', 'options', s2)
[c2, gof2] = fit (T_2_results_x', T_2_results_y', h2, 'problem', t_delay)

plot(c2, 'k');
xlabel ('t [sec]');
ylabel ('M_{trans} [Volt]');
text ('Position', [3e-3, 0.4], 'String', {strcat(' Fit results: M_0=', num2str(c2.M0)), strcat(' {T_2}=', num2str(c2.T_2))});

% for a comparison we'll try to perform fit for the linear function after
% logarithm evauation
figure;
plot (T_2_results_x, log(T_2_results_y), 'o');
title ('T_2 calculation out of 8% CP experiment (linear fitting curve)')
legend ('experimental data')
hold on

%s4 = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0, 0],'Upper', [Inf, Inf], 'Startpoint',[1e-6, max_Voltage]);
h5 = fittype ('poly1');
[c5 gof5] = fit (T_2_results_x', log(T_2_results_y)', h5)
plot(c5, 'k');
xlabel ('t [sec]');
ylabel ('ln(M_{trans}) [Volt]');
text ('Position', [2e-3, -0.8], 'String', {strcat(' Fit results: M_0=', num2str(exp(c5.p2))), strcat(' {T_2}=', num2str(-1/c5.p1))});




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% T_2 Calculation out of CPMG experiment%%%%%%

[C I] = max (T2_8percent_CPMG(:,2));
max_Voltage = C;
t_delay = T2_8percent_CPMG(I, 1);
figure;
plot(T2_8percent_CPMG(:,1), T2_8percent_CPMG(:, 2));
deriv_Volt_mg = diff(T2_8percent_CPMG(:, 2))./diff(T2_8percent_CPMG(:,1)); % we obtained the differential
deriv_Volt_mg = [0; deriv_Volt_mg];
figure;
plot(T2_8percent_CPMG(:,1), deriv_Volt);

deriv_threshed_mg = (deriv_Volt_mg>Threshold);
d_deriv_border_mg = find(diff (deriv_threshed_mg)>0)
j = 1;
for i=1:2:length(d_deriv_border_mg)-1
    
    [C I] = max (T2_8percent_CPMG(d_deriv_border_mg(i):d_deriv_border_mg(i+1), 2));
    T_2_results_mg_y(j) = C;
    T_2_results_mg_x(j) = T2_8percent_CP(I+d_deriv_border_mg(i), 1);
    j = j+1;
end

figure;
plot (T_2_results_mg_x, T_2_results_mg_y);
title ('T_2 calculation out of 8% CPMG experiment');
legend ('CPMG experimental data');
hold on

s3 = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0, 0],'Upper', [Inf, Inf], 'Startpoint',[1e-6, max_Voltage]);
h3 = fittype ('M0*exp(-(t-t_del)/T_2)', 'coefficients', {'T_2', 'M0'}, 'independent', 't', 'problem', 't_del', 'options', s3)
[c3, gof3] = fit (T_2_results_mg_x', T_2_results_mg_y', h3, 'problem', t_delay)

plot(c3, 'k');
xlabel ('t [sec]');
ylabel ('M_{trans} [Volt]');
text ('Position', [3e-3, 0.4], 'String', {strcat(' Fit results: M_0=', num2str(c3.M0)), strcat(' {T_2}=', num2str(c3.T_2))});


% for a comparison we'll try to perform fit for the linear function after
% logarithm evauation
figure;
plot (T_2_results_mg_x, log(T_2_results_mg_y), 'o');
title ('T_2 calculation out of 8% CPMG experiment (linear fitting curve)')
legend ('experimental data')
hold on

%s4 = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0, 0],'Upper', [Inf, Inf], 'Startpoint',[1e-6, max_Voltage]);
h4 = fittype ('poly1');
[c4 gof4] = fit (T_2_results_mg_x', log(T_2_results_mg_y)', h4)
plot(c4, 'k');
xlabel ('t [sec]');
ylabel ('ln(M_{trans}) [Volt]');
text ('Position', [2e-3, -0.8], 'String', {strcat(' Fit results: M_0=', num2str(exp(c4.p2))), strcat(' {T_2}=', num2str(-1/c4.p1))});

