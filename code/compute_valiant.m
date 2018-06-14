function [t] = compute_valiant(null_dist,epsilon,counts)
d = length(null_dist);
n = sum(counts);
null_cdf = cumsum(null_dist);
tail_min = find(null_cdf >= (1 - epsilon/8),1);
tail = (tail_min + 1):d;

% bulk indices
bulk = setxor(2:d,tail);

t(1) = sum(counts(tail));
t(2) = 0;
for j = 1:length(bulk)
    t(2) = t(2) + ((counts(bulk(j)) - n*null_dist(bulk(j)))^2 - counts(bulk(j)))/(null_dist(bulk(j))^(2/3));
end
end