%% 选择复制
%选择
function [newpop]=selection(pop,fitvalue)
totalfit=sum(fitvalue);%求适应值之和
fitvalue=fitvalue/totalfit;%单个个体被选择的概率
[px,py]=size(pop);
fitvalue=cumsum(fitvalue);
ms=sort(rand(px,1));%从小到大排列
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