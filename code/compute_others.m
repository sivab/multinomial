function [t] = compute_others(null_dist,counts)


% 1: Chi-squared statistic
% 2: L2 statistic
% 3: L1 statistic
% 4: Wasserman statistic
% 5: LRT
d = length(null_dist);
n = sum(counts);
t = zeros(1,5);
    for i = 1:d
        t(1) = t(1) + ((counts(i) - n*null_dist(i))^2 - n*null_dist(i))/null_dist(i);
        t(2) = t(2) + ((counts(i) - n*null_dist(i))^2 - n*null_dist(i));
        t(3) = t(3) + abs(counts(i) - n*null_dist(i));
        t(4) = t(4) + ((counts(i) - n*null_dist(i))^2 - counts(i))/(max(1,d*null_dist(i)));
        if (counts(i) > 0)
            t(5) = t(5) + (counts(i)*log(counts(i)/(n*null_dist(i))))/n;
        end
    end 
end

