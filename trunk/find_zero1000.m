function [beg] = find_zero1000(A_array)

for i=1:100:length(A_array)-1000
    if (max(abs(A_array(i:i+1000, 2)))<0.05)
        beg =i;
        return;
    end
end
    beg = -1;
end
