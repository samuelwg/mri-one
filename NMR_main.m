%%% NMR Experiment main program %%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

%%% Gyromagnetic Ratio calculation

B0 = [3.6 3.56 3.59 3.60 3.53];
freq = [(15.27257+15.27198)/2 (15.29292+15.29306)/2 (15.29660+15.29764)/2 (15.30097+15.30002)/2 (15.30392+15.30260)/2]
gamma = freq./B0;
conc = [8 4 1 0.5 0.25];
figure;
plot(conc, gamma, 'o')
title('Gyromagnetic Ratio as function of concentration');
xlabel('Concentration, perc');
ylabel('\gamma [mHz/kG]');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T1_list = dir ('T1*.csv');
for i=1:length(T1_list)
    [T1_results(i,1), T1_results(i,2), T1_results(i,3), T1_results(i,4)] = T1(T1_list(i).name);
end

%%% we are going to plot T_1 as function of concentration
figure;
plot (T1_results(:,4), T1_results(:,2));
title ('T_1 as function of concentration of CuSo_4 ');
xlabel('Concentration, perc');
ylabel('T_1, [msec]');


T2_star_list = dir ('T2_star*.csv');

for i=1:length(T2_star_list)
    [T2_star_results(i,1), T2_star_results(i,2), T2_star_results(i,3), T2_star_results(i,4)] = T2_star(T2_star_list(i).name);
end

%%% we are going to plot T_1 as function of concentration
figure;
plot (T2_star_results(:,4), T2_star_results(:,2), 'o');
title ('{T_2}^* as function of concentration of CuSo_4 ');
xlabel('Concentration, perc');
ylabel('{T_2}^*, [sec]');




T2_list = dir ('TT2_*CPMG.csv');

for i=1:length(T2_list)
    [T2_results_CPMG(i,1), T2_results_CPMG(i,2), T2_results_CPMG(i,3), T2_results_CPMG(i,4)] = T2(T2_list(i).name, 1);
end

%% Now we'll plot T2 as a function of concentration
figure;
plot (T2_results_CPMG(:,4), T2_results_CPMG(:,2), '-o');
title ('{T_2} (CPMG sequence) as function of concentration of CuSo_4 ');
xlabel('Concentration, perc');
ylabel('{T_2}, [sec]');

T2_CP_list = dir ('TT2_*CP.csv');

for i=1:length(T2_CP_list)
    [T2_results_CP(i,1), T2_results_CP(i,2), T2_results_CP(i,3), T2_results_CP(i,4)] = T2(T2_CP_list(i).name, 0);
end
%% Now we'll plot T2 as a function of concentration
figure;
plot (T2_results_CP(:,4), T2_results_CP(:,2), '-o');
title ('{T_2} (CP sequence) as function of concentration of CuSo_4 ');
xlabel('Concentration, perc');
ylabel('{T_2}, [sec]');