function [G_k,G_best_route,G_best_length,K,best_route,best_length,length_ave]=VRP_original(C,NC_max,m,Alpha,Beta,Rho,Q,W,t) %#ok<INUSD,STOUT>
%%
%--------------------------------------------------------------------------
%%参数说明
%%G_K:各代最优车辆数目
%%G_best_route：各代最优路线
%%G_best_length：各代最优路线长度
%%K：最优车辆数目
%%best_route：最优路线
%%best_length：最优路线长度
%%length_ave：各代平均长度
%%C：DC和工厂的坐标
%%NC_max：最大迭代次数
%%m：蚂蚁数目
%%Alpha：重要度系数
%%Beta：能见度系数
%%Rho：挥发度系数
%%Q：信息更新参数
%%W:车辆载重量
%%
%%第一步 初始化变量和系数

m=10;Alpha=2;Beta=1;gama=1;Rho=0.1;NC_max=100;Q=15;W=9;qq=0.05; % #ok<NASGU>
C=[0 0
    0 -1
    0 3
    -2 -2
    -3 -3
    3 -1
    -4 0
    -4 -1
    1 -2
    1 -1
    1 3
    3 4
    -3 0
    2 0
    1 -3
    2 -1
    2 1
    1 -4
    -3 2
    -1 -1]; %  #ok<NASGU>
t=[0 1.5 1.8 2 0.8 1.5 1.0 2.5 3.0 1.7 0.6 0.2 2.4 1.9 2.0 0.7 0.5 2.2 3.1 0.1];

%构造仓库和工厂等之间的距离矩阵

n=size(C,1);%n表示问题的规模（城市个数） 
D=zeros(n,n);%D表示完全图的赋权邻接矩阵 
for i=1:n 
for j=1:n 
if i~=j 
D(i,j)=((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5; 
else 
D(i,j)=eps; 
end 
D(j,i)=D(i,j); 
end 
end 

%构造节省量矩阵
U=zeros(n,n);%U表示工厂之间的连接和其与仓库之间连接能节省的距离
for i=1:n 
for j=1:n 
if i~=j 
U(i,j)=D(i,1)+D(j,1)-D(i,j);
else 
U(i,j)=eps; 
end 
U(j,i)=U(i,j); 
end 
end 
load_w=0;
Eta=1./D;%Eta为启发因子，这里设为距离的倒数 
Tau=ones(n,n);%Tau为信息素矩阵 
Tabu=zeros(m,n);%存储并记录路径的生成 
NC=1;%迭代计数器 
G_best_route=[NC_max,n];%各代最佳路线 
G_best_length=inf.*ones(NC_max,1);%各代最佳路线的长度 
length_ave=zeros(NC_max,1);%各代路线的平均长度

%%第二步，把蚂蚁放到DC内
while NC<=NC_max%停止条件之一：达到最大迭代次数 
    Tabu(:,1)=randi([1,1],m,1);
%%第三步，m只蚂蚁按照要求的方法选择工厂，并完成周游 
for i=1:m
    visited=Tabu(i,:);
    visited=visited(visited>0);
    to_visit=setdiff(1:n,visited);
    c_temp=length(to_visit);
    j=1;
       while j<=n
          if ~isempty(to_visit)
    %visit_thisant=[];
    %while to_visit~=[1]
    %for ee=1:(2*c_temp) 
   %while c_temp>=2
   %if length(to_visit)>1
    %visited=Tabu(i,:);
    %to_visit=setdiff(1:n,visited);
    %to_visit=[1,to_visit];
    %c_temp=length(to_visit);
    %visit_thisant=[];
   % if c_temp~=1
        %for b=1:c_temp
            %if (load_w+t(to_visit(b)))<=W
                %visit_thisant=[visit_thisant,to_visit(b)];
            %end
        %end
        %jj=visited(visited>0);
%% 按照规则选下一个工厂或者是回到仓库
%for k=1:length(visit_thisant)
    for k=1:length(to_visit)
  x(k)=(Tau(visited(end),to_visit(k))^Alpha)*(Eta(visited(end),to_visit(k))^Beta)*(U(visited(end),to_visit(k))^gama);
    end   
%         ww=rand;
% if ww<qq
%     Select=find(max(x));
%     %Tabu(i,length(visited(visited>0))+1)=to_visit(Select(1)); 
% else
x=x/(sum(x)); 
%按概率原则选取下一个城市 
xcum=cumsum(x); 
Select=find(xcum>=rand);
%Tabu(i,length(visited(visited>0))+1)=to_visi(Select(1)); 
% end
if isempty(Select)
    Select=1;
    load_w=load_w+t(Select);
else
load_w=load_w+t(to_visit(Select(1)));
end
if load_w>W
    Select=1;
       j=j-1;
    load_w=0;
    Tabu(i,length(visited)+1)=Select(1);
else
Tabu(i,length(visited)+1)=to_visit(Select(1)); 
end
%if Select(1)==1

%else
    %load_w=load_w+t(Select(1));
%end
          end
    visited=Tabu(i,:);
    visited=visited(visited>0);
    to_visit=setdiff(1:n,visited);
    x=[];
            if visited(end)~=1
   Tabu(i,1:(length(visited)+1))=[visited,1];
            end
        j=j+1;
       end
    load_w=0;
end
%% 第四步记录本代各种参数
L=zeros(m,1); 
for i=1:m 
MM=Tabu(i,:); 
R=MM(MM>0);
for j=1:(length(R)-1)
L(i)=L(i)+D(R(j),R(j+1)); 
end 
end 

G_best_length(NC)=min(L); 
pos=find(L==G_best_length(NC)); 
G_best_route(NC,1:length(Tabu(pos(1),:)))=Tabu(pos(1),:);
%%应用2-opt方法对最优解进行更新
select=find(G_best_route(NC,:)==1);
FF=[];
LM=0;
for a=1:(length(select)-1)
   y_G_best_route=G_best_route(NC,select(a):select(a+1));
   al=length(y_G_best_route);
   T=zeros((length(select)-1),1);
   for d=1:(al-1)
       T(a)=T(a)+D(y_G_best_route(d),y_G_best_route(d+1));
   end
   for b=2:(al-1)
       for c=(b+1):(al-1)
           DD=y_G_best_route;
           temp1=DD(b);
           temp2=DD(c);
           DD(b)=temp2;
           DD(c)=temp1;
           TT=zeros(1);
           for d=1:(al-1)
       TT=TT+D(DD(d),DD(d+1));
           end
           if TT<T(a)
               T(a)=TT;
            y_G_best_route=DD;
           end
       end
   end
   if a>=2
       y_G_best_route=y_G_best_route(2:al);
   end
   FF=[FF,y_G_best_route];
   LM=LM+T(a);
end
   G_best_length(NC)=LM;
   G_best_route(NC,1:length(FF))=FF;
   FF=[];
   LM=0;
   %%2-opt优化完毕
length_ave(NC)=mean(L); 
NC=NC+1;
%% 第五步更新信息素
Delta_Tau=zeros(n,n); 
for i=1:m 
MM=Tabu(i,:); 
R=MM(MM>0);
for j=1:(length(R)-1) 
Delta_Tau(R(j),R(j+1))=Delta_Tau(R(j),R(j+1))+Q/L(i); 
end 
end 
Tau=(1-Rho).*Tau+Delta_Tau;
%% 第六步：禁忌表清零 
Tabu=zeros(m,n); 
load_w=0;
end
%%第七步：输出结果 
Pos=find(G_best_length==min(G_best_length)); 
best_route=G_best_route(Pos(1),:);
best_length=G_best_length(Pos(1));
best_route=best_route(best_route>0);
best_route
best_length
%% 第八步：绘制散点图和巡游过程图

figure(1)
plot([C(best_route,1)],[C(best_route,2)],'-*')
 figure(2)
plot(G_best_length) 
hold on 
plot(length_ave) 

end
