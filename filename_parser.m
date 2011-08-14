function [ sol_perc_num ] = filename_parser( file_name )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%%%% Parsing the filename %%%%%%%%%
[mat idx] = regexp(file_name, '\d+', 'match', 'start');
sol_perc = mat(2);
sol_perc_num=str2num(char(sol_perc));
switch sol_perc_num
    case 5
        sol_perc_num = 0.5;
    case 25
        sol_perc_num =0.25;
end

end

