function [ threshold_val] = find_thresh( deriv_array )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[C, I] = max (deriv_array);
deriv_array(I) = 0;
[C, I] = max (deriv_array);
deriv_array(I) = 0;
[C, I] = max (deriv_array);
threshold_val = 0.4*C; 
end

