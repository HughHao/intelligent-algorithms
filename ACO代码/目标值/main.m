clear all
clc
clf
popsize=20;
chromlength=10;
pc=0.6;
pm=0.001;

pop=fu(popsize,chromlength);
for i=1:20
    [objvalue]=calobjvalue(pop);%计算目标函数
    fitvalue=calfitvalue(objvalue);
    [newpop]=selection(pop,fitvalue);
    [newpop]=crossover(pop,pc);
    [newpop]=mutation(pop,pm);
    [bestindividual,bestfit]=best(pop,fitvalue);
    y(i)=max(bestfit);
    n(i)=i;
    pop5=bestindividual;
    x(i)=decodechrom(pop5,1,chromlength)*10/1023;
    pop=newpop;
end
fplot('10*sin(5*x)+7*cos(4*x)',[0,10])
hold on 
plot(x,y,'r*')
hold off
[z index]=max(y);
%计算最大西对应的x
x5=x(index)
y=z