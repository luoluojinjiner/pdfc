
tic;  % 开始计时
clc
clear
close all
BX=xlsread('F:\笔记本备份\F\The fivth paper()\双参数收敛情况\Jfuntion2.xlsx');
h=plot(BX);
[ymin,k] = min(BX);

set(gcf,'PaperType','a3');
 set(h,'LineWidth',1.7);
xlabel('Step','FontSize',18,'Fontname','Times newman');
ylabel('J','FontSize',18,'Fontname','Times newman');
set(gca,'FontName','Times New Roman','FontSize',18);
axis([2,250,min(BX)-1,max(BX)+1]);

