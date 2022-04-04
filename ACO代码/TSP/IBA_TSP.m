function IBA_TSP
 
CityNum=31; %you can choose 10, 30, 50, 75
[dislist,Clist]=tsp(CityNum);%距离清单，城市列单
 
inn=31; %初始种群大小
gnmax=500; %最大代数
pc=0.8; %交叉概率
pm=0.8; %变异概率
 
%产生初始种群
zhongqun=zeros(inn,CityNum);
for i=1:inn
    zhongqun(i,:)=randperm(CityNum);%城市数量大小的随机乱序序列
end
[~,p]=objf(zhongqun,dislist);%目标适应度值
 
gn=1;
ymean=zeros(gn,1);
ymax=zeros(gn,1);
xmax=zeros(inn,CityNum);
scnew=zeros(inn,CityNum);
smnew=zeros(inn,CityNum);
while gn<gnmax+1
   for j=1:2:inn
      seln=sel(p);  %选择操作
      scro=cro(zhongqun,seln,pc);  %交叉操作
      scnew(j,:)=scro(1,:);
      scnew(j+1,:)=scro(2,:);
      smnew(j,:)=mut(scnew(j,:),pm);  %变异操作
      smnew(j+1,:)=mut(scnew(j+1,:),pm);
   end
   zhongqun=smnew;  %产生了新的种群
   [f,p]=objf(zhongqun,dislist);  %计算新种群的适应度
   %记录当前代最好和平均的适应度
   [fmax,nmax]=max(f);
   ymean(gn)=1000/mean(f);
   ymax(gn)=1000/fmax;
   %记录当前代的最佳个体
   x=zhongqun(nmax,:);
   xmax(gn,:)=x;
   drawTSP(Clist,x,ymax(gn),gn,0);
   gn=gn+1;
end
[min_ymax,jn]=min(ymax);
drawTSP(Clist,xmax(jn,:),min_ymax,jn,1);
 
figure(2);
plot(ymax,'r'); hold on;
plot(ymean,'b');grid on;
title('搜索过程');
legend('最优解','平均解');
fprintf('遗传算法得到的最短距离:%.2f\n',min_ymax);
fprintf('遗传算法得到的最短路线');
disp(xmax(jn,:));
end
 
%------------------------------------------------
%计算所有种群的适应度
function [f,p]=objf(s,dislist)
 
inn=size(s,1);  %读取种群大小
f=zeros(inn,1);
for i=1:inn
   f(i)=CalDist(dislist,s(i,:));  %计算函数值，即适应度
end
f=1000./f'; %取距离倒数
%根据个体的适应度计算其被选择的概率
fsum=0;
for i=1:inn
   fsum=fsum+f(i)^15;% 让适应度越好的个体被选择概率越高
end
ps=zeros(inn,1);
for i=1:inn
   ps(i)=f(i)^15/fsum;
end
 
%计算累积概率
p=zeros(inn,1);
p(1)=ps(1);
for i=2:inn
   p(i)=p(i-1)+ps(i);
end
p=p';
end
 
%--------------------------------------------------
%根据变异概率判断是否变异
function pcc=pro(pc)
test(1:100)=0;
l=round(100*pc);
test(1:l)=1;
n=round(rand*99)+1;
pcc=test(n);   
end
 
%--------------------------------------------------
%“选择”操作
function seln=sel(p)
 
seln=zeros(2,1);
%从种群中选择两个个体，最好不要两次选择同一个个体
for i=1:2
   r=rand;  %产生一个随机数
   prand=p-r;
   j=1;
   while prand(j)<0
       j=j+1;
   end
   seln(i)=j; %选中个体的序号
   if i==2&&j==seln(i-1)    %%若相同就再选一次
       r=rand;  %产生一个随机数
       prand=p-r;
       j=1;
       while prand(j)<0
           j=j+1;
       end
   end
end
end
 
%------------------------------------------------
%“交叉”操作
function scro=cro(s,seln,pc)
 
bn=size(s,2);
pcc=pro(pc);  %根据交叉概率决定是否进行交叉操作，1则是，0则否
scro(1,:)=s(seln(1),:);
scro(2,:)=s(seln(2),:);
if pcc==1
   c1=round(rand*(bn-2))+1;  %在[1,bn-1]范围内随机产生一个交叉位
   c2=round(rand*(bn-2))+1;
   chb1=min(c1,c2);
   chb2=max(c1,c2);
   middle=scro(1,chb1+1:chb2);
   scro(1,chb1+1:chb2)=scro(2,chb1+1:chb2);
   scro(2,chb1+1:chb2)=middle;
   for i=1:chb1 %似乎有问题
       while find(scro(1,chb1+1:chb2)==scro(1,i))
           zhi=find(scro(1,chb1+1:chb2)==scro(1,i));
           y=scro(2,chb1+zhi);
           scro(1,i)=y;
       end
       while find(scro(2,chb1+1:chb2)==scro(2,i))
           zhi=find(scro(2,chb1+1:chb2)==scro(2,i));
           y=scro(1,chb1+zhi);
           scro(2,i)=y;
       end
   end
   for i=chb2+1:bn
       while find(scro(1,1:chb2)==scro(1,i))
           zhi=logical(scro(1,1:chb2)==scro(1,i));
           y=scro(2,zhi);
           scro(1,i)=y;
       end
       while find(scro(2,1:chb2)==scro(2,i))
           zhi=logical(scro(2,1:chb2)==scro(2,i));
           y=scro(1,zhi);
           scro(2,i)=y;
       end
   end
end
end
 
%--------------------------------------------------
%“变异”操作
function snnew=mut(snew,pm)
 
bn=size(snew,2);
snnew=snew;
 
pmm=pro(pm);  %根据变异概率决定是否进行变异操作，1则是，0则否
if pmm==1
   c1=round(rand*(bn-2))+1;  %在[1,bn-1]范围内随机产生一个变异位
   c2=round(rand*(bn-2))+1;
   chb1=min(c1,c2);
   chb2=max(c1,c2);
   x=snew(chb1+1:chb2);
   snnew(chb1+1:chb2)=fliplr(x);
end
end
 
%------------------------------------------------
%城市位置坐标
function [dislist,cityn]=tsp(n)
dislist=zeros(n,n);
if n==10
    city10=[0.4 0.4439;0.2439 0.1463;0.1707 0.2293;0.2293 0.761;0.5171 0.9414;
        0.8732 0.6536;0.6878 0.5219;0.8488 0.3609;0.6683 0.2536;0.6195 0.2634];%10 cities d'=2.691
    for i=1:10
        for j=1:10
            dislist(i,j)=((city10(i,1)-city10(j,1))^2+(city10(i,2)-city10(j,2))^2)^0.5;
        end
    end
    cityn=city10;
end
if n==31
    city30=[1304.0	2312.0;3639.0	1315.0;4177.0	2244.0;3712.0	1399.0;3488.0	1535.0;3326.0	1556.0;3238.0	1229.0;
4196.0	1004.0;4312.0	790.0;4386.0	570.0;3007.0	1970.0;2562.0	1756.0;2788.0	1491.0;2381.0	1676.0;1332.0	695.0;
3715.0	1678.0;3918.0	2179.0;4061.0	2370.0;3780.0	2212.0;3676.0	2578.0;4029.0	2838.0;4263.0	2931.0;
3429.0	1908.0
3507.0	2367.0
3394.0	2643.0
3439.0	3201.0
2935.0	3240.0
3140.0	3550.0
2545.0	2357.0
2778.0	2826.0
2370.0	2975.0
];%30 cities d'=423.741 by D B Fogel
    for i=1:31
        for j=1:31
            dislist(i,j)=((city30(i,1)-city30(j,1))^2+(city30(i,2)-city30(j,2))^2)^0.5;
        end
    end
    cityn=city30;
end
 
if n==50
    city50=[31 32;32 39;40 30;37 69;27 68;37 52;38 46;31 62;30 48;21 47;25 55;16 57;
        17 63;42 41;17 33;25 32;5 64;8 52;12 42;7 38;5 25; 10 77;45 35;42 57;32 22;
        27 23;56 37;52 41;49 49;58 48;57 58;39 10;46 10;59 15;51 21;48 28;52 33;
        58 27;61 33;62 63;20 26;5 6;13 13;21 10;30 15;36 16;62 42;63 69;52 64;43 67];%50 cities d'=427.855 by D B Fogel
    for i=1:50
        for j=1:50
            dislist(i,j)=((city50(i,1)-city50(j,1))^2+(city50(i,2)-city50(j,2))^2)^0.5;
        end
    end
    cityn=city50;
end
 
if n==75
    city75=[48 21;52 26;55 50;50 50;41 46;51 42;55 45;38 33;33 34;45 35;40 37;50 30;
        55 34;54 38;26 13;15 5;21 48;29 39;33 44;15 19;16 19;12 17;50 40;22 53;21 36;
        20 30;26 29;40 20;36 26;62 48;67 41;62 35;65 27;62 24;55 20;35 51;30 50;
        45 42;21 45;36 6;6 25;11 28;26 59;30 60;22 22;27 24;30 20;35 16;54 10;50 15;
        44 13;35 60;40 60;40 66;31 76;47 66;50 70;57 72;55 65;2 38;7 43;9 56;15 56;
        10 70;17 64;55 57;62 57;70 64;64 4;59 5;50 4;60 15;66 14;66 8;43 26];%75 cities d'=549.18 by D B Fogel
    for i=1:75
        for j=1:75
            dislist(i,j)=((city75(i,1)-city75(j,1))^2+(city75(i,2)-city75(j,2))^2)^0.5;
        end
    end
    cityn=city75;
end
end
 
%------------------------------------------------
%适应度函数
function F=CalDist(dislist,s)
 
DistanV=0;
n=size(s,2);
for i=1:(n-1)
    DistanV=DistanV+dislist(s(i),s(i+1));
end
DistanV=DistanV+dislist(s(n),s(1));
F=DistanV;
 
end
 
%------------------------------------------------
%画图
function drawTSP(Clist,x,y,p,index)
CityNum=size(Clist,1);
for i=1:CityNum-1
    plot([Clist(x(i),1),Clist(x(i+1),1)],[Clist(x(i),2),Clist(x(i+1),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
    text(Clist(x(i),1),Clist(x(i),2),['  ',int2str(x(i))]);
    text(Clist(x(i+1),1),Clist(x(i+1),2),['  ',int2str(x(i+1))]);
    hold on;
end
plot([Clist(x(CityNum),1),Clist(x(1),1)],[Clist(x(CityNum),2),Clist(x(1),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
title([num2str(CityNum),'城市TSP']);
if index==0&&CityNum~=10
    text(1000,1000,['第 ',int2str(p),' 代','  最短距离为 ',num2str(y)]);
else
    text(1000,1000,['最终搜索结果：最短距离 ',num2str(y),'， 在第 ',num2str(p),' 代达到']);
end
if CityNum==10
    if index==0
        text(1000,1000,['第 ',int2str(p),' 代','  最短距离为 ',num2str(y)]);
    else
        text(1000,1000,['最终搜索结果：最短距离 ',num2str(y),'， 在第 ',num2str(p),' 代达到']);
    end
end
hold off;
pause(0.001); 
end
