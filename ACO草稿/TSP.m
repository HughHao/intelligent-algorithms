%% I. 清空环境变量
clear all
clc

%% II. 导入数据
load citys.mat

%% III. 计算城市间相互距离
n = size(citys,1);
D = 1000*ones(n,n);
D(13,19) = sqrt(sum((citys(13,:) - citys(19,:)).^2));
D(19,13) = D(13,19);
for i = 1:4
    for j = 1:4
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = 5:8
    for j = 5:8
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = 9:13
    for j = 9:13
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = 14:19
    for j = 14:19
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
D(20,18) = sqrt(sum((citys(20,:) - citys(18,:)).^2));
D(18,20) = D(20,18);
for i = 21:24
    for j = 21:24
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = 25:27
    for j = 1:4
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = 28:31
    for j = 28:31
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = 32:34
    for j = 32:34
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = 35:36
    for j = 35:36
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = 37:38
    for j = 37:38
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = 39:40
    for j = 39:40
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = [1 9 14 21 28]
    for j = [1 9 14 21 28]
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = [15 22 26 29 33 36 39]
    for j = [15 22 26 29 33 36 39]
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = [2 6 10 16 23 27 30 34 37 40]
    for j = [2 6 10 16 23 27 30 34 37 40]
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = [3 7 11 17 24 31 38]
    for j = [3 7 11 17 24 31 38]
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
for i = [4 8 12]
    for j = [4 8 12]
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end
%% IV. 初始化参数
m = 40;                              % 蚂蚁数量
alpha = 1;                           % 信息素重要程度因子
beta = 5;                            % 启发函数重要程度因子
rho = 0.1;                           % 信息素挥发因子
Q = 1;                               % 常系数
Eta = 1./D;                          % 启发函数
Tau = ones(n,n);                     % 信息素矩阵
Table = zeros(m,n);                  % 路径记录表
iter = 1;                            % 迭代次数初值
iter_max = 200;                      % 最大迭代次数 
Route_best = zeros(iter_max,n);      % 各代最佳路径       
Length_best = zeros(iter_max,1);     % 各代最佳路径的长度  
Length_ave = zeros(iter_max,1);      % 各代路径的平均长度  
tabu = 1000*ones(n,n);
for i = 1:n
    for j = 1:n
        if D(i,j)>=1000
            tabu(i,j) = j;
        end
    end
end
%% V. 迭代寻找最佳路径
while iter <= iter_max
     % 随机产生各个蚂蚁的起点城市
      start = zeros(m,1);
      for i = 1:m
          %temp = randperm(n);%随机打乱一个数字序列，其内的参数决定了随机数的范围。
          start(i) = 20;
      end
      Table(:,1) = start;%将start赋值到Table的第一列 
      citys_index = 1:n;
      % 逐个蚂蚁路径选择
      jl=1;
      for i = 1:m
          % 逐个城市路径选择
         while Table(i,jl)~=20
             jl=jl+1;
             %tabu = Table(i,1:(j - 1)); % 已访问的城市集合(禁忌表),Table的第i行，1到j-1列
             allow_index = ~ismember(citys_index,tabu(i,:));%结果是0，1矩阵
             allow = citys_index(allow_index);  % 待访问的城市集合
             P = allow;%重置P，每只蚂蚁都有新P
             % 计算城市间转移概率
             for k = 1:length(allow)%允许访问的城市数量
                 P(k) = Tau(Table(i,jl-1),allow(k))^alpha * Eta(Table(i,jl - 1),allow(k))^beta;
             end%上面Tau是1矩阵，Eta是1/D的β次方
             P = P/sum(P);%sum(A)，若A为行向量时，不指定dim或指定dim为2，则自动计算成所有行向量数值的和，如果指定dim为1，则计算结果为一个行向量，且与原来的行向量相同。
             % 轮盘赌法选择下一个访问城市，轮盘赌法即
             Pc = cumsum(P); %cumsum作用是累加P第一列概率，最后一个为1    
            target_index = find(Pc >= rand); %找出Pc中所有大于等于rand的代号
            target = allow(target_index(1));%第一个目标是allow中target首位
            Table(i,jl) = target;%第i行第j列为新目标
         end
      end
      % 计算各个蚂蚁的路径距离
      Len = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          for j = 1:(length(Route) - 1)
              Len(i) = Len(i) + D(Route(j),Route(j + 1));
          end
          Len(i) = Len(i) + D(Route(end),Route(1));
      end
      % 计算最短路径距离及平均距离
      if iter == 1
          [min_Length,min_index] = min(Len);%求出最短路径的值和其在length中的位置
          Length_best(iter) = min_Length;  
          Length_ave(iter) = mean(Len);
          Route_best(iter,:) = Table(min_index,:);%求出最短路径为table中的最小值代码行
      else
          [min_Length,min_index] = min(Len);
          Length_best(iter) = min(Length_best(iter - 1),min_Length);
          Length_ave(iter) = mean(Len);
          if Length_best(iter) == min_Length
              Route_best(iter,:) = Table(min_index,:);
          else
              Route_best(iter,:) = Route_best((iter-1),:);
          end
      end
      % 更新信息素
      Delta_Tau = zeros(n,n);%初始矩阵ΔTau=零
      % 逐个蚂蚁计算
      for i = 1:m
          % 逐个城市计算
          for j = 1:(length(Route) - 1)%控制循环到倒数第二个城市
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Len(i);
          end%从现在城市到下个城市的ΔTau=Q/length(i)加上上一个ΔTau
          Delta_Tau(Table(i,end),Table(i,1)) = Delta_Tau(Table(i,end),Table(i,1)) + Q/Len(i);
      end%该行最后一个城市到第一个城市的ΔTau另外公式
      Tau = (1-rho) * Tau + Delta_Tau;%Tau的现在信息素浓度，减去衰减部分
    % 迭代次数加1，清空路径记录表
    iter = iter + 1;
    %Table = zeros(m,n);
end

%% VI. 结果显示
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['最短距离:' num2str(Shortest_Length)]);
disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);%先给出最短路线的代号，再连接上第一个城市代号。中间没有分号隔开

%% VII. 绘图
figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...%citys(Shortest_Route,1)是所有城市的横坐标，然后citys(Shortest_Route(1),1)]是第一个城市的横坐标，
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'ko-');%citys(Shortest_Route,2);citys(Shortest_Route(1),2)是所有城市的纵坐标与第一个城市的纵坐标。与plot([x,y])作用类似
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       起点');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       终点');
xlabel('城市位置横坐标')
ylabel('城市位置纵坐标')
title(['蚁群算法优化路径(最短距离:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r:')
legend('最短距离','平均距离')
xlabel('迭代次数')
ylabel('距离')
title('蚁群算法各代最短距离与平均距离对比')