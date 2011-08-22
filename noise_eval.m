function [offset] = noise_eval(file_name)
noise_sig = load (file_name);
noise_beg = length(noise_sig)-2000;
noise_sig_part = noise_sig(noise_beg:end, 2);
%figure;
%subplot(4,1,1);
% plot(noise_sig_part, '.');

[c, x] = hist (noise_sig_part)
% subplot(4,1,2);
% bar(x,c);

[c1, i1] = max(c);
max_noise_val = x(i1)
offset = max_noise_val;
end