%% ѡ����
%ѡ��
function [newpop]=selection(pop,fitvalue)
totalfit=sum(fitvalue);%����Ӧֵ֮��
fitvalue=fitvalue/totalfit;%�������屻ѡ��ĸ���
[px,py]=size(pop);
fitvalue=cumsum(fitvalue);
ms=sort(rand(px,1));%��С��������
fitin=1;
newin=1;
while newin<=px
    if (ms(newin))<fitvalue(fitin)
        newpop(newin)=pop(fitin);
        newin=newin+1;
    else
        fitin=fitin+1;
    end
end