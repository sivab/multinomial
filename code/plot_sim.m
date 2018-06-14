function [] = plot_sim(data,labels,fname,num_repeats)

% Format of data
% first row corresponds to x points
% remaining rows are corresponding y values
% Format of labels: 
% row 1: title
% row 2: xlabel
% row 3: ylabel
% row 4 onwards legend

nlines = size(data,1) - 1;

if (nlines >= 6)
    fprintf('Things are going to crash!');
end

c = {'-r+',':b*','-.kd','--go','-yx','--cs'};

cmtemp = colormap(lines(6));
cm = zeros(size(cmtemp));
cm(1,:) = cmtemp(1,:);
cm(2,:) = cmtemp(6,:);% n-> number of lines to plot
cm(3,:) = cmtemp(2,:);
cm(4,:) = cmtemp(5,:);

close all;
figure;
%set(gca,'fontsize',18);
%set(gca,'fontname','Timesnewroman');
set(gca, ... 
  'FontSize', 16, ... 
  'XColor'      , [.1 .1 .1], ... 
  'YColor'      , [.1 .1 .1], ... 
  'LineWidth'   , 2);

hold on;

title(labels{1},'Interpreter','latex');
grid on;
%set(gca,'GridLineStyle','--');
hold on;


i = 1;
errorbar(data(1,:),data(i+1,:),sqrt(data(i+1,:)/num_repeats),c{i},'LineWidth',2,'MarkerSize',10,'Color',cm(i,:));
i = 4;
errorbar(data(1,:),data(i+1,:),sqrt(data(i+1,:)/num_repeats),c{i},'LineWidth',2,'MarkerSize',10,'Color',cm(i,:));
i = 2;
errorbar(data(1,:),data(i+1,:),sqrt(data(i+1,:)/num_repeats),c{i},'LineWidth',2,'MarkerSize',10,'Color',cm(i,:));
i = 3;
errorbar(data(1,:),data(i+1,:),sqrt(data(i+1,:)/num_repeats),c{i},'LineWidth',2,'MarkerSize',10,'Color',cm(i,:));

xlim([0 max(data(1,:))]);
ylim([0 1]);

labelstemp = labels;
labels{5} = labelstemp{7};
labels{6} = labelstemp{5};
labels{7} = labelstemp{6};
xlabel(labels{2},'Interpreter','latex');
ylabel(labels{3},'Interpreter','latex');
legend(labels{4:end},'Location','best');
print(fname,'-depsc2');

end


