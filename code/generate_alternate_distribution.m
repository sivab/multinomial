function [alt_dist] = generate_alternate_distribution(sim_num,null_dist,eps_distance)

% Sim num:
% 1: perturb every coordinate uniformly
% 2: perturb only the first two coordinates
% 3: perturb every coordinate proportional to +/- p_i^{2/3}
% 4: perturb every coordinate proportional to +/- p_i
d = length(null_dist);
alt_dist = zeros(d,1);

v = ones(d,1);
t = randperm(d);
q = floor(d/2);
v(t(1:q)) = -1;
if (2*q ~= d)
    v(t(d),1) = 0;
end

if (sim_num == 1)
    % alt_dist = null_dist + v*eps_distance/(2*q);
    % doing it Larry's way
    for i = 1:d
        if (mod(i,2) == 1)
            alt_dist(i) = null_dist(i) + eps_distance/d;
        else
            alt_dist(i) = null_dist(i) - eps_distance/d;
        end
    end
elseif (sim_num == 2)
    %s = 2;
    alt_dist = null_dist;
    alt_dist(1) = null_dist(1) + eps_distance/2;
    alt_dist(2) = null_dist(2) + eps_distance/2;
    %     for i = 1:s
    %         alt_dist(i) = null_dist(i) + sign(rand-0.5)*eps_distance/s;
    %     end
elseif (sim_num == 3)
    pert_vec = null_dist.^(2/3);
    pert_vec = eps_distance*pert_vec/sum(abs(pert_vec));
    pert_vec = min(pert_vec,null_dist);
    alt_dist = null_dist + v.*pert_vec;
elseif (sim_num == 4)
    pert_vec = null_dist;
    pert_vec = eps_distance*pert_vec/sum(abs(pert_vec));
    alt_dist = null_dist + v.*pert_vec;
end

% cleanup

alt_dist(alt_dist < 0) = 0;
alt_dist = alt_dist/sum(alt_dist);


end