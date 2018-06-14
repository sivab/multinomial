function [tchisq] = gof_tchisq(X,null_dist)
d = length(X);
n = sum(X);
tchisq = 0;
for i = 1:d
    tchisq = tchisq + ((X(i) - n*null_dist(i))^2 - X(i))/max(1,d*null_dist(i));
end
end