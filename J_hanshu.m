
tic;  % ��ʼ��ʱ
clc
clear
close all
BX=xlsread('F:\�ʼǱ�����\F\The fivth paper()\˫�����������\Jfuntion2.xlsx');
h=plot(BX);
[ymin,k] = min(BX);

set(gcf,'PaperType','a3');
 set(h,'LineWidth',1.7);
xlabel('Step','FontSize',18,'Fontname','Times newman');
ylabel('J','FontSize',18,'Fontname','Times newman');
set(gca,'FontName','Times New Roman','FontSize',18);
axis([2,250,min(BX)-1,max(BX)+1]);

