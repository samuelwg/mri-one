%%% NMR Experiment main program %%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
 pos = get(0, 'ScreenSize');

%%% Gyromagnetic Ratio calculation

B0 = [3.6 3.56 3.59 3.60 3.53];
freq = [(15.27257+15.27198)/2 (15.29292+15.29306)/2 (15.29660+15.29764)/2 (15.30097+15.30002)/2 (15.30392+15.30260)/2]
gamma = freq./B0;
conc = [8 4 1 0.5 0.25];
h1 = figure('Color', 'white');
plot(conc, gamma, 'o')
title('Gyromagnetic Ratio as function of concentration');
xlabel('Concentration, perc');
ylabel('\gamma [mHz/kG]');
export_fig(h1, 'gmagnet.pdf');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T1_list = dir ('T1*.csv');
h0=figure('PaperOrientation', 'landscape','position', pos, 'Color', 'white');
for i=1:length(T1_list)
    [T1_results(i,1), T1_results(i,2), T1_results(i,3), T1_results(i,4)] = T1(T1_list(i).name, i);
end

%%% we are going to plot T_1 as function of concentration
subplot(2,3, i+1);
plot (T1_results(:,4), T1_results(:,2));
title ('T_1 as function of concentration of CuSo_4 ');
xlabel('Concentration, perc');
ylabel('T_1, [msec]');
export_fig(h0, 'T1res.pdf')

T2_star_list = dir ('T2_star*.csv');
h2 = figure('PaperOrientation', 'landscape','position', pos, 'Color', 'white');
for i=1:length(T2_star_list)
    [T2_star_results(i,1), T2_star_results(i,2), T2_star_results(i,3), T2_star_results(i,4)] = T2_star(T2_star_list(i).name, i, 0);
end

%%% we are going to plot T_2_star as function of concentration
subplot (2,3,i+1)
plot (T2_star_results(:,4), T2_star_results(:,2), 'o');
title ('{T_2}^* as function of concentration of CuSo_4 ');
xlabel('Concentration, perc');
ylabel('{T_2}^*, [sec]');
export_fig (h2, 'T2starres.pdf');



T2_list = dir ('TT2_*CPMG.csv');
h3 = figure('PaperOrientation', 'landscape','position', pos, 'Color', 'white');
for i=1:length(T2_list)
    [T2_results_CPMG(i,1), T2_results_CPMG(i,2), T2_results_CPMG(i,3), T2_results_CPMG(i,4)] = T2(T2_list(i).name, 1, i, 0);
end

%% Now we'll plot T2 as a function of concentration
subplot (2,3, i+1);
plot (T2_results_CPMG(:,4), T2_results_CPMG(:,2), '-o');
title ('{T_2} (CPMG sequence) as function of concentration of CuSo_4 ');
xlabel('Concentration, perc');
ylabel('{T_2}, [sec]');
export_fig(h3, 'T2_CPMGres.pdf')

T2_CP_list = dir ('TT2_*CP.csv');
h4 = figure('PaperOrientation', 'landscape','position', pos, 'Color', 'white');
for i=1:length(T2_CP_list)
    [T2_results_CP(i,1), T2_results_CP(i,2), T2_results_CP(i,3), T2_results_CP(i,4)] = T2(T2_CP_list(i).name, 0, i, 0);
end
%% Now we'll plot T2 as a function of concentration
subplot (2,3, i+1);
plot (T2_results_CP(:,4), T2_results_CP(:,2), '-o');
title ('{T_2} (CP sequence) as function of concentration of CuSo_4 ');
xlabel('Concentration, perc');
ylabel('{T_2}, [sec]');
export_fig(h4,'T2_CPres.pdf')

Ts_list = dir ('*T2*.csv');
for i = 1:length(Ts_list)
    a(i) = noise_eval(Ts_list(i).name)
end


%%% now we will compare results with an moving average (3 & 5 points) filter %%%%%

T2_list = dir ('TT2_*CPMG.csv');
h7 = figure('PaperOrientation', 'landscape','position', pos, 'Color', 'white');
for i=1:length(T2_list)
    [T2_results_CPMG_f3(i,1), T2_results_CPMG_f3(i,2), T2_results_CPMG_f3(i,3), T2_results_CPMG_f3(i,4)] = T2_filt(T2_list(i).name, 1, i, 3);
end

%% Now we'll plot T2 as a function of concentration
subplot (2,3, i+1);
plot (T2_results_CPMG(:,4), T2_results_CPMG(:,2), '-o');
title ('{T_2} (CPMG sequence) as function of concentration of CuSo_4 fitered 3p');
xlabel('Concentration, perc');
ylabel('{T_2}, [sec]');
export_fig(h7, 'T2_CPMGres_filt3.pdf')



T2_list = dir ('TT2_*CPMG.csv');
h8 = figure('PaperOrientation', 'landscape','position', pos, 'Color', 'white');
for i=1:length(T2_list)
    [T2_results_CPMG_f5(i,1), T2_results_CPMG_f5(i,2), T2_results_CPMG_f5(i,3), T2_results_CPMG_f5(i,4)] = T2_filt(T2_list(i).name, 1, i, 5);
end

%% Now we'll plot T2 as a function of concentration
subplot (2,3, i+1);
plot (T2_results_CPMG(:,4), T2_results_CPMG(:,2), '-o');
title ('{T_2} (CPMG sequence) as function of concentration of CuSo_4 fitered 5p');
xlabel('Concentration, perc');
ylabel('{T_2}, [sec]');
export_fig(h8, 'T2_CPMGres_filt5.pdf')

h9 = figure('Color', 'white');
plot(T2_results_CPMG(:,4), T2_results_CPMG(:,2), '-or');
hold on;
plot(T2_results_CPMG_f3(:,4), T2_results_CPMG_f3(:,2), '-ok');
plot(T2_results_CPMG_f5(:,4), T2_results_CPMG_f5(:,2), '-ob');
legend('not filtered data', 'moving average 3 points', 'moving average 5 points')
title ('{T_2} (CPMG sequence) as function of concentration of CuSo_4 with various filter length');
xlabel('Concentration, perc');
ylabel('{T_2}, [sec]');
export_fig(h9, 'T2_CPMGres_filt_diff.pdf')

h10 = figure('Color', 'white');

plot(T2_results_CPMG(:,4), T2_results_CPMG(:,3), '-or');
hold on;
plot(T2_results_CPMG_f3(:,4), T2_results_CPMG_f3(:,3), '-ok');
plot(T2_results_CPMG_f5(:,4), T2_results_CPMG_f5(:,3), '-ob');
legend('not filtered data', 'moving average 3 points', 'moving average 5 points', 'Location', 'SE')
title ('{T_2} (CPMG sequence) fit goodness as function of concentration of CuSo_4 with various filter length');
xlabel('Concentration, perc');
ylabel('R^2');
export_fig(h10, 'T2_CPMGres_filt_diff_Rsq.pdf')



T2_star_list = dir ('T2_star*.csv');
h11 = figure('PaperOrientation', 'landscape','position', pos, 'Color', 'white');
for i=1:length(T2_star_list)
    [T2_star_results(i,1), T2_star_results(i,2), T2_star_results(i,3), T2_star_results(i,4)] = T2_star(T2_star_list(i).name, i, 7);
end

%%% we are going to plot T_2_star as function of concentration
subplot (2,3,i+1)
plot (T2_star_results(:,4), T2_star_results(:,2), 'o');
title ('{T_2}^* as function of concentration of CuSo_4 filtered 3p');
xlabel('Concentration, perc');
ylabel('{T_2}^*, [sec]');
export_fig (h11, 'T2starres_f3.pdf');

