function [chisq] = gof_chisq(X,null_dist)
d = length(X);
n = sum(X);
chisq = 0;
for i = 1:d
    if (null_dist(i) > 0)
     chisq = chisq + ((X(i) - n*null_dist(i))^2 - n*null_dist(i))/(n*null_dist(i));
    end
end

end