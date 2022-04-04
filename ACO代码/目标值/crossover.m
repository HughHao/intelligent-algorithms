%% ½»²æ
function [newpop]=crossover(pop,pc)
[px,py]=size(pop);
newpop=ones(size(pop));
for i =1:2:px-1
    if (rand<pc)
        cpoint=round(rand*py);
        newpop(i,:)=[pop(i,1:cpoint),pop(i+1,cpoint+1:py)];
        newpop(i,:)=[pop(i,1:cpoint),pop(i,cpoint+1:py)];
    else
        newpop(i,:)=pop(i);
        newpop(i,:)=pop(i);
    end
end