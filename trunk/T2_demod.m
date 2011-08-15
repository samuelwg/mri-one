function [ M0, T2_res, R_square percentage ] = T2_demod( file_name, CPMG )

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
t0_point_ind = t0_point(1);
T2_demod = demod(T2(t0_point_ind :end, 2), 833.3e2 ,2.5e6, 'qam' );
figure;
plot (T2_demod);
