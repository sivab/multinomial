function [] = paper1()

% This file plots the null distribution of goodness-of-fit 
% test statistics

n = 400;
d = 1000;
p = zeros(1,d);
N = 1000;

for i = 1:d
    p(i) = 1/i^2;
end
p = p/sum(p);

% null distribution
stat_t = zeros(N,1);
stat_c = zeros(N,1);
for i = 1:N
    X = mnrnd(n,p);
    stat_t(i) = gof_tchisq(X,p);
    stat_c(i) = gof_chisq(X,p);
end

stat_t = sort(stat_t,'ascend');

close all
figure
histogram(stat_t,50,'Normalization','probability','EdgeColor','w','FaceColor',[0 0.5 0.5]);
set(gca,'fontsize',18);
set(gca,'fontname','Timesnewroman');
hold on;
grid on;
set(gca, ... 
  'FontSize', 16, ... 
  'XColor'      , [.1 .1 .1], ... 
  'YColor'      , [.1 .1 .1], ... 
  'LineWidth'   , 2);
title('Histogram of Truncated $\chi^2$ Statistic','Interpreter','latex');
print('Plots/hist_trunc.eps','-depsc2');


figure
stat_t = stat_t/std(stat_t);
qqplot(stat_t);
set(gca,'fontsize',18);
set(gca,'fontname','Timesnewroman');
hold on;
grid on;
set(gca, ... 
  'FontSize', 16, ... 
  'XColor'      , [.1 .1 .1], ... 
  'YColor'      , [.1 .1 .1], ... 
  'LineWidth'   , 2);
title('QQ Plot of Truncated $\chi^2$ Statistic','Interpreter','latex');
xlabel('Standard Normal Quantiles','Interpreter','latex');
ylabel('Sample Quantiles','Interpreter','latex');
print('Plots/qq_tchisq.eps','-depsc2');

figure
histogram(stat_c,50,'Normalization','probability','EdgeColor','w','FaceColor',[0 0.5 0.5]);
set(gca,'fontsize',18);
set(gca,'fontname','Timesnewroman');
hold on;
grid on;
set(gca, ... 
  'FontSize', 16, ... 
  'XColor'      , [.1 .1 .1], ... 
  'YColor'      , [.1 .1 .1], ... 
  'LineWidth'   , 2);
title('Histogram of Null Distribution of Classical $\chi^2$ Statistic','Interpreter','latex');
print('Plots/hist_chisq.eps','-depsc2');


figure
stat_c = stat_c/std(stat_c);
qqplot(stat_c);
set(gca,'fontsize',18);
set(gca,'fontname','Timesnewroman');
hold on;
grid on;
set(gca, ... 
  'FontSize', 16, ... 
  'XColor'      , [.1 .1 .1], ... 
  'YColor'      , [.1 .1 .1], ... 
  'LineWidth'   , 2);
title('QQ Plot for Classical $\chi^2$ Statistic','Interpreter','latex');
xlabel('Standard Normal Quantiles','Interpreter','latex');
ylabel('Sample Quantiles','Interpreter','latex');
print('Plots/qq_chisq.eps','-depsc2');     
end

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

function [tchisq] = gof_tchisq(X,null_dist)
d = length(X);
n = sum(X);
tchisq = 0;
for i = 1:d
    tchisq = tchisq + ((X(i) - n*null_dist(i))^2 - X(i))/max(1,d*null_dist(i));
end
end