%% I. ��ջ�������
clear all
clc

%% II. ��������
load citys_data.mat

%% III. ������м��໥����
n = size(citys,1);
D = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;      
        end
    end    
end

%% IV. ��ʼ������
m = 35;                              % ��������
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

%% V. ����Ѱ�����·��
while iter <= iter_max
     % ��������������ϵ�������
      start = zeros(m,1);
      for i = 1:m
          temp = randperm(n);%�������һ���������У����ڵĲ���������������ķ�Χ��
          start(i) = temp(1);
      end
      Table(:,1) = start;%��start��ֵ��Table�ĵ�һ�� 
      citys_index = 1:n;
      % �������·��ѡ��
      for i = 1:m
          % �������·��ѡ��
         for j = 2:n
             tabu = Table(i,1:(j - 1)); % �ѷ��ʵĳ��м���(���ɱ�),Table�ĵ�i�У�1��j-1��
             allow_index = ~ismember(citys_index,tabu);%�����0��1����
             allow = citys_index(allow_index);  % �����ʵĳ��м���
             P = allow;%����P��ÿֻ���϶�����P
             % ������м�ת�Ƹ���
             for k = 1:length(allow)%������ʵĳ�������
                 P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
             end%����Tau��1����Eta��1/D�Ħ´η�
             P = P/sum(P);%sum(A)����AΪ������ʱ����ָ��dim��ָ��dimΪ2�����Զ������������������ֵ�ĺͣ����ָ��dimΪ1���������Ϊһ��������������ԭ������������ͬ��
             % ���̶ķ�ѡ����һ�����ʳ��У����̶ķ���
             Pc = cumsum(P); %cumsum�������ۼ�P��һ�и��ʣ����һ��Ϊ1    
            target_index = find(Pc >= rand); %�ҳ�Pc�����д��ڵ���rand�Ĵ���
            target = allow(target_index(1));%��һ��Ŀ����allow��target��λ
            Table(i,j) = target;%��i�е�j��Ϊ��Ŀ��
         end
      end
      % ����������ϵ�·������
      Length = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);
          for j = 1:(n - 1)
              Length(i) = Length(i) + D(Route(j),Route(j + 1));
          end
          Length(i) = Length(i) + D(Route(n),Route(1));
      end
      % �������·�����뼰ƽ������
      if iter == 1
          [min_Length,min_index] = min(Length);%������·����ֵ������length�е�λ��
          Length_best(iter) = min_Length;  
          Length_ave(iter) = mean(Length);
          Route_best(iter,:) = Table(min_index,:);%������·��Ϊtable�е���Сֵ������
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
      % ������Ϣ��
      Delta_Tau = zeros(n,n);%��ʼ����Tau=��
      % ������ϼ���
      for i = 1:m
          % ������м���
          for j = 1:(n - 1)%����ѭ���������ڶ�������
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
          end%�����ڳ��е��¸����еĦ�Tau=Q/length(i)������һ����Tau
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
      end%�������һ�����е���һ�����еĦ�Tau���⹫ʽ
      Tau = (1-rho) * Tau + Delta_Tau;%Tau��������Ϣ��Ũ�ȣ���ȥ˥������
    % ����������1�����·����¼��
    iter = iter + 1;
    Table = zeros(m,n);
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