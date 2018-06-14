function [] = paper5()

% unif null sparse perturbation

alt_distance = 0:0.05:1;
n = 400;
d = 1000;
num_repeats = 1000;

null_dist = ones(d,1);
null_dist = null_dist/sum(null_dist);



x = zeros(length(alt_distance),1);
y = zeros(length(alt_distance),6);
err = zeros(length(alt_distance),6);


for i = 1:length(alt_distance)
    succ_v = zeros(num_repeats,1);
    succ_o = zeros(num_repeats,5);
    alt_dist = generate_alternate_distribution(2,null_dist,alt_distance(i));
    act_distance = sum(abs(alt_dist - null_dist));
    if (i == 1)
        [tv,to] = null_thresholds(null_dist,n,act_distance);
    else
        [tv,~] = null_thresholds(null_dist,n,act_distance);
    end
    
    for j = 1:num_repeats
        counts = mnrnd(n,alt_dist);
        v = compute_valiant(null_dist,act_distance,counts);
        t = compute_others(null_dist,counts);
        
        if ((v(1) > tv(1))||(v(2) > tv(2)))
            succ_v(j) = 1;
        end
        
        for k = 1:5
            if (t(k) >= to(k))
                succ_o(j,k) = 1;
            end
        end              
    end
    x(i) = act_distance;
    y(i,1:5) = sum(succ_o,1)/num_repeats;
    err(i,1:5) = (sum(succ_o,1)/num_repeats).*(1 - (sum(succ_o,1)/num_repeats));
    y(i,6) = sum(succ_v)/num_repeats;
    err(i,6) = y(i,6)*(1-y(i,6));
    
end

c = {'-r+','--go',':b*','-.kd','-yx','--cs'};

cm = colormap(lines(6)); % n-> number of lines to plot

close all;
figure;
hold on;



set(gca, ... 
  'FontSize', 16, ... 
  'XColor'      , [.1 .1 .1], ... 
  'YColor'      , [.1 .1 .1], ... 
  'LineWidth'   , 2);

cmind = [1,5,6,2,3,4];
for i = 1:6
    %plot(x,y(:,i),c{i},'LineWidth',3,'MarkerSize', 12);
    errorbar(x,y(:,i),sqrt(err(:,i)/num_repeats),c{i},'LineWidth',2,'MarkerSize',10,'Color',cm(cmind(i),:));
end



hLegend = legend('Chi-sq.','L2','L1','Truncated chi-sq','LRT','2/3-rd and tail','Location','east');
set(hLegend  ,'FontSize'   , 16);
grid on;

title('Uniform Null, Sparse Alternate','Interpreter','latex');
xlim([0 max(x)]);
ylim([0 1]);
xlabel('$$\ell_1$$ Distance','Interpreter','latex');
ylabel('Power','Interpreter','latex');
set(gcf, 'PaperPositionMode', 'auto');
filename = sprintf('Plots/n%d_d%d_unif_sparse_full.eps',n,d);
print(filename,'-depsc2');

end