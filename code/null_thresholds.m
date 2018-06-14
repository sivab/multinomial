function [valiant_threshold, other_thresholds,temp_tail,temp_bulk] = null_thresholds(null_dist,n,epsilon)

N = 1000;
alpha = 0.05;

tv = zeros(N,2);
to = zeros(N,5);

for i = 1:N
    counts = mnrnd(n,null_dist);
    tv(i,:) = compute_valiant(null_dist,epsilon,counts);
    to(i,:) = compute_others(null_dist,counts);
end

other_thresholds = zeros(5,1);

for i = 1:5
    temp = sort(to(:,i),'ascend');
    other_thresholds(i) = temp(N*(1-alpha) - 1);
end

temp_tail = sort(tv(:,1),'ascend');
temp_bulk = sort(tv(:,2),'ascend');

valiant_threshold(1) = temp_tail(N*(1-alpha/2) - 1);
valiant_threshold(2) = temp_bulk(N*(1-alpha/2) - 1);



end
