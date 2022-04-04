%test_lianxu_function
%第一层和最后一层只有0
%中间层都是0-9  len
%T蚂蚁的路径矩阵，起始点是第一层的0，终点是第d+2层的0
%tao全部d+2层的信息素矩阵，代表每一层上每一节点的经过概率
len=10;
tao0=0.01;
d=7;
tao=ones(d+2,len,len);
tao=tao*tao0;%所有蚂蚁共用
rou=0.8;
N0=20;
T=zeros(N0,d+2);%N0只蚂蚁，d+2层，赋予初始值为0，也就是说所有蚂蚁第一层次的都从第一个城市出发
Q0=0.8;
d=7;
iteration=100;
N0=20;
alpha=0.8;

max_result=0;
max_route=[];
max_x=0;
for i=1:iteration%迭代次数，迭代1000次，也就是所有蚂蚁走过1000遍
    %a=1;%第一层从数组下标1开始
    for j=2:d+1%层数，从初始点第一层0出发，到最后一层的0结束
        %b=1;
        for k=1:N0%蚂蚁编码
            b=select_k_city(tao,Q0,k,j,T,d);
            %b=select_k_city(tao,Q0,j,a,d);%求出第k只蚂蚁，第j层节点位置
            T(k,j)=b(1)-1;%T数组保存的是第k蚂蚁第j层的具体数值，应该为城市编号减1
            
            %执行局部更新城市之间的信息素
            %k    第k层
            %tao  两层之间的信息素
            %rou  比例
            %tao0 剩余信息素0.01
            %T    每一只蚂蚁下一步到达的城市.第一维是第几只蚂蚁、第二位是目前是第几步
            %n    第几只蚂蚁
            tao=update_tao(tao,rou,tao0,T,j,k);
        end
        %a=b;
    end
    
    %当所有蚂蚁全局都跑了一遍后，进行全局性更新
    %7) 根据公式(4) ～ (6) 评选出最优蚂蚁并执行全局更新规则;
    %共有n只蚂蚁
    %tao层与层之间的信息素
    %d是小数点后精确到d位

    %X所有蚂蚁的具体数值数组
    %idx最小的蚂蚁编码，也就是数组下标
    %alpha是一常量
    %[tao,X,idx]=update_all_tao(tao,n,alpha,d)
    [tao,X,idx]=update_all_tao(tao,N0,alpha,d,T);
    
    
    max_cur=myfunction(X(idx(1)));
    
    if max_cur>max_result
       max_result=max_cur;
       max_route=T(idx,:);
       max_x=X(idx);        
    end
    
    
end
max_result
max_route
max_x
