%% I. ��ջ�������
clear all
clc

%% II. ��������
load citys.mat

%% III. ������м��໥����
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
%% IV. ��ʼ������
m = 40;                              % ��������
alpha = 1;                           % ��Ϣ����Ҫ�̶�����
beta = 5;                            % ����������Ҫ�̶�����
rho = 0.1;                           % ��Ϣ�ػӷ�����
Q = 1;                               % ��ϵ��
Eta = 1./D;                          % ��������
Tau = ones(n,n);                     % ��Ϣ�ؾ���
Table = zeros(m,n);                  % ·����¼��
iter = 1;                            % ����������ֵ
iter_max = 200;                      % ���������� 
Route_best = zeros(iter_max,n);      % �������·��       
Length_best = zeros(iter_max,1);     % �������·���ĳ���  
Length_ave = zeros(iter_max,1);      % ����·����ƽ������  
tabu = 1000*ones(n,n);
for i = 1:n
    for j = 1:n
        if D(i,j)>=1000
            tabu(i,j) = j;
        end
    end
end
%% V. ����Ѱ�����·��
while iter <= iter_max
     % ��������������ϵ�������
      start = zeros(m,1);
      for i = 1:m
          %temp = randperm(n);%�������һ���������У����ڵĲ���������������ķ�Χ��
          start(i) = 20;
      end
      Table(:,1) = start;%��start��ֵ��Table�ĵ�һ�� 
      citys_index = 1:n;
      % �������·��ѡ��
      jl=1;
      for i = 1:m
          % �������·��ѡ��
         while Table(i,jl)~=20
             jl=jl+1;
             %tabu = Table(i,1:(j - 1)); % �ѷ��ʵĳ��м���(���ɱ�),Table�ĵ�i�У�1��j-1��
             allow_index = ~ismember(citys_index,tabu(i,:));%�����0��1����
             allow = citys_index(allow_index);  % �����ʵĳ��м���
             P = allow;%����P��ÿֻ���϶�����P
             % ������м�ת�Ƹ���
             for k = 1:length(allow)%������ʵĳ�������
                 P(k) = Tau(Table(i,jl-1),allow(k))^alpha * Eta(Table(i,jl - 1),allow(k))^beta;
             end%����Tau��1����Eta��1/D�Ħ´η�
             P = P/sum(P);%sum(A)����AΪ������ʱ����ָ��dim��ָ��dimΪ2�����Զ������������������ֵ�ĺͣ����ָ��dimΪ1���������Ϊһ��������������ԭ������������ͬ��
             % ���̶ķ�ѡ����һ�����ʳ��У����̶ķ���
             Pc = cumsum(P); %cumsum�������ۼ�P��һ�и��ʣ����һ��Ϊ1    
            target_index = find(Pc >= rand); %�ҳ�Pc�����д��ڵ���rand�Ĵ���
            target = allow(target_index(1));%��һ��Ŀ����allow��target��λ
            Table(i,jl) = target;%��i�е�j��Ϊ��Ŀ��
         end
      end
      % ����������ϵ�·������
      Len = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          for j = 1:(length(Route) - 1)
              Len(i) = Len(i) + D(Route(j),Route(j + 1));
          end
          Len(i) = Len(i) + D(Route(end),Route(1));
      end
      % �������·�����뼰ƽ������
      if iter == 1
          [min_Length,min_index] = min(Len);%������·����ֵ������length�е�λ��
          Length_best(iter) = min_Length;  
          Length_ave(iter) = mean(Len);
          Route_best(iter,:) = Table(min_index,:);%������·��Ϊtable�е���Сֵ������
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
      % ������Ϣ��
      Delta_Tau = zeros(n,n);%��ʼ����Tau=��
      % ������ϼ���
      for i = 1:m
          % ������м���
          for j = 1:(length(Route) - 1)%����ѭ���������ڶ�������
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Len(i);
          end%�����ڳ��е��¸����еĦ�Tau=Q/length(i)������һ����Tau
          Delta_Tau(Table(i,end),Table(i,1)) = Delta_Tau(Table(i,end),Table(i,1)) + Q/Len(i);
      end%�������һ�����е���һ�����еĦ�Tau���⹫ʽ
      Tau = (1-rho) * Tau + Delta_Tau;%Tau��������Ϣ��Ũ�ȣ���ȥ˥������
    % ����������1�����·����¼��
    iter = iter + 1;
    %Table = zeros(m,n);
end

%% VI. �����ʾ
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['��̾���:' num2str(Shortest_Length)]);
disp(['���·��:' num2str([Shortest_Route Shortest_Route(1)])]);%�ȸ������·�ߵĴ��ţ��������ϵ�һ�����д��š��м�û�зֺŸ���

%% VII. ��ͼ
figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...%citys(Shortest_Route,1)�����г��еĺ����꣬Ȼ��citys(Shortest_Route(1),1)]�ǵ�һ�����еĺ����꣬
     [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'ko-');%citys(Shortest_Route,2);citys(Shortest_Route(1),2)�����г��е����������һ�����е������ꡣ��plot([x,y])��������
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       ���');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       �յ�');
xlabel('����λ�ú�����')
ylabel('����λ��������')
title(['��Ⱥ�㷨�Ż�·��(��̾���:' num2str(Shortest_Length) ')'])
figure(2)
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r:')
legend('��̾���','ƽ������')
xlabel('��������')
ylabel('����')
title('��Ⱥ�㷨������̾�����ƽ������Ա�')