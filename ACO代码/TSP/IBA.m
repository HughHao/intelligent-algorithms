%% I. 清空环境变量
clear all
clc

%% II. 导入数据
load citys_data.mat

%% III. 计算城市间相互距离
n = size(citys,1);
D = zeros(n,n);
for i = 1:n
    for iter = 1:n
        if i ~= iter
            D(i,iter) = sqrt(sum((citys(i,:) - citys(iter,:)).^2));
        else
            D(i,iter) = 1e-4;      
        end
    end    
end

%% 初始化参数
m=40;
Eta = 1./D;      
iter_max=500;
Table = zeros(m,n);                  % 路径记录表
Route_best = zeros(iter_max,n);      % 各代最佳路径       
Length_best = zeros(iter_max,1);     % 各代最佳路径的长度  
Length_ave = zeros(iter_max,1);      % 各代路径的平均长度  
for iter=1:iter_max
% 随机产生各个种群的起点城市
      start = zeros(m,1);
     for i = 1:m
          temp = randperm(n);%随机打乱一个数字序列，其内的参数决定了随机数的范围。
          start(i) = temp(1);
      end
      Table(:,1) = start;%将start赋值到Table的第一列 
      citys_index = 1:n;
      % 逐个蚂蚁路径选择
      for i = 1:m
          % 逐个城市路径选择
         for j = 2:n
             tabu = Table(i,1:(j - 1)); % 已访问的城市集合(禁忌表),Table的第i行，1到j-1列
             allow_index = ~ismember(citys_index,tabu);%结果是0，1矩阵
             allow = citys_index(allow_index);  % 待访问的城市集合
             P = allow;%重置P，每只蚂蚁都有新P
             % 计算城市间转移概率
             for k = 1:length(allow)%允许访问的城市数量
                 P(k) =  Eta(tabu(end),allow(k));
             end%上面Tau是1矩阵，Eta是1/D的β次方
             P = P/sum(P);%sum(A)，若A为行向量时，不指定dim或指定dim为2，则自动计算成所有行向量数值的和，如果指定dim为1，则计算结果为一个行向量，且与原来的行向量相同。
             % 轮盘赌法选择下一个访问城市，轮盘赌法即
             Pc = cumsum(P); %cumsum作用是累加P第一列概率，最后一个为1    
            target_index = find(Pc >= rand); %找出Pc中所有大于等于rand的代号
            target = allow(target_index(1));%第一个目标是allow中target首位
            Table(i,j) = target;%第i行第j列为新目标
         end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % 计算各个蚂蚁的路径距离
      Length = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          for j = 1:(n - 1)
              Length(i) = Length(i) + D(Route(j),Route(j + 1));
          end
          Length(i) = Length(i) + D(Route(n),Route(1));
      end
      % 计算最短路径距离及平均距离
      if iter == 1
          [min_Length,min_index] = min(Length);%求出最短路径的值和其在length中的位置
          Length_best(iter) = min_Length;  
          Length_ave(iter) = mean(Length);
          Route_best(iter,:) = Table(min_index,:);%求出最短路径为table中的最小值代码行
      else
          [min_Length,min_index] = min(Length);
          Length_best(iter) = min(Length_best(iter - 1),min_Length);
          Length_ave(iter) = mean(Length);
          if Length_best(iter) == min_Length
              Route_best(iter,:) = Table(min_index,:);
          else
              Route_best(iter,:) = Route_best((iter-1),:);
          end
      end
end
%% VI. 结果显示
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['最短距离:' num2str(Shortest_Length)]);
disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);%先给出最短路线的代号，再连接上第一个城市代号。中间没有分号隔开

%% VII. 绘图
figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...%citys(Shortest_Route,1)是所有城市的横坐标，然后citys(Shortest_Route(1),1)]是第一个城市的横坐标，
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'o-');%citys(Shortest_Route,2);citys(Shortest_Route(1),2)是所有城市的纵坐标与第一个城市的纵坐标。与plot([x,y])作用类似
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       起点');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       终点');
xlabel('城市位置横坐标')
ylabel('城市位置纵坐标')
title(['牵制平衡算法优化路径(最短距离:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r:')
legend('最短距离','平均距离')
xlabel('迭代次数')
ylabel('距离')
title('各代最短距离与平均距离对比')