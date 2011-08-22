function [offset] = noise_eval(file_name)
noise_sig = load (file_name);
zero_beg = find_zero1000(noise_sig)
%noise_beg = length(noise_sig)-1000;
noise_sig_part = noise_sig(zero_beg:zero_beg+1000, 2);
%figure;
%subplot(4,1,1);
% plot(noise_sig_part, '.');

[c, x] = hist (noise_sig_part);
%figure;
%subplot(4,1,2);
%bar(x,c);

[c1, i1] = max(c);
max_noise_val = x(i1)
offset = max_noise_val;
end