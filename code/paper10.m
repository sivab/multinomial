function [] = paper10()

% we simulate P power law, Q is a sparse perturbation

n_1 = 200;
n_2 = 200;
d = 400;
N = 1000;
eps = 0:0.05:1;

p = zeros(1,d);
for i = 1:d
    p(i) = 1/i;
end
p = p/sum(p);
save_data = zeros(length(eps),3);

fchisq = @ts_chisq;
fellone = @ts_ellone;
ftchisq = @gof_tchisq;
felltwo = @ts_elltwo;

for outer = 1:length(eps)
    outer
    q = p;
    q(1) = p(1) + eps(outer)/2;
    q(2) = p(2) + eps(outer)/2;
    q = q/sum(q);
    
    succ_chisq = 0;
    succ_ellone = 0;
    succ_tchisq = 0;
    succ_elltwo = 0;
    
    gof_cutoff = sim_cutoff(p,n_2,ftchisq);
    for i = 1:N
        
        X = mnrnd(n_1,p);
        Y = mnrnd(n_2,q);
        
        chisq_cut = perm_cutoff(X,Y,fchisq);
        ellone_cut = perm_cutoff(X,Y,fellone);
        elltwo_cut = perm_cutoff(X,Y,felltwo);
        
        chisq_act = ts_chisq(X,Y);
        ellone_act = ts_ellone(X,Y);
        elltwo_act = ts_elltwo(X,Y);
        tchisq_act = gof_tchisq(Y,p);
        
        
        if (chisq_act >= chisq_cut)
            succ_chisq = succ_chisq + (1/N);
        end
        if (ellone_act >= ellone_cut)
            succ_ellone = succ_ellone + (1/N);
        end
        if (tchisq_act >= gof_cutoff)
            succ_tchisq = succ_tchisq + (1/N);
        end
        if (elltwo_act >= elltwo_cut)
            succ_elltwo = succ_elltwo + (1/N);
        end
    end
    save_data(outer,1) = succ_chisq;
    save_data(outer,2) = succ_ellone;
    save_data(outer,3) = succ_tchisq;
    save_data(outer,4) = succ_elltwo;
end
%fprintf('\n Chi Sq.: %d Ell One: %d, GoF: %d, distance: %d',succ_chisq,succ_ellone,succ_tchisq,sum(abs(p-q)));

data(1,:) = eps;
data(2:5,:) = save_data';
labels{1} = 'Two-Sample, Power Law';
labels{2} = '$\ell_1$ Distance';
labels{3} = 'Power';
labels{4} = 'Chi-Sq.';
labels{5} = 'L1';
labels{6} = 'Oracle GoF';
labels{7} = 'L2';
plot_sim(data,labels,'Plots/two_sample_power_balance.eps',N);

end

function [cutoff] = sim_cutoff(null_dist,sample_size,test)

B = 1000;
bootstat = zeros(B,1);

for i = 1:B
    S = mnrnd(sample_size,null_dist);
    bootstat(i) = test(S,null_dist);
end

bootstat = sort(bootstat,'ascend');
cutoff = bootstat(0.95*B);

end

function [cutoff] = perm_cutoff(X,Y,test)

B = 500;
expX = [];
expY = [];
d = length(X);
for i = 1:d
    expX = [expX i*ones(1,X(i))];
    expY = [expY i*ones(1,Y(i))];
end
expXY = [expX expY];
n_1 = size(expX,2);
n = size(expXY,2);

bootstat = zeros(B,1);

newX = zeros(1,d);
newY = zeros(1,d);

for i = 1:B
    arr = expXY(randperm(n));
    newexpX = arr(1,1:n_1);
    newexpY = arr(1,(n_1 + 1):end);
    newX = accumarray(newexpX',1);
    newX = [newX;zeros(d-size(newX,1),1)];
    
    newY = accumarray(newexpY',1);
    newY = [newY;zeros(d-size(newY,1),1)];
    
    %     for j = 1:d
    %         newX(j) = sum(newexpX == j);
    %         newY(j) = sum(newexpY == j);
    %     end
    bootstat(i) = test(newX,newY);
end

bootstat = sort(bootstat,'ascend');
cutoff = bootstat(B*0.95);

end

function [tchisq] = gof_tchisq(X,null_dist)
d = length(X);
n = sum(X);
tchisq = 0;
for i = 1:d
    tchisq = tchisq + ((X(i) - n*null_dist(i))^2 - X(i))/max(1,d*null_dist(i));
end
end

function [chisq] = ts_chisq(X,Y)

n_1 = sum(X);
n_2 = sum(Y);
temp = ((n_2*X - n_1*Y).^2 - (n_2^2)*X - (n_1^2)*Y)./(X + Y);
temp(isnan(temp)) = 0;
chisq = sum(temp);

end

function [ellone] = ts_ellone(X,Y)

n_1 = sum(X);
n_2 = sum(Y);

ellone = sum(abs(X/n_1 - Y/n_2));
end

function [elltwo] = ts_elltwo(X,Y)

n_1 = sum(X);
n_2 = sum(Y);

elltwo = sum((X/n_1 - Y/n_2).^2);
end


