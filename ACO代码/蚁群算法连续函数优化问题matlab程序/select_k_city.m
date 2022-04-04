%{
5) 根据公式(1) 和(2) 选择蚂蚁n在第k 层应该到达的城市;
tao  任意两个城市之间的信息素，三个维度，第一维表示第几层，第二维表示
Q0   是一个常数，主要是作为轮盘赌法的临界点
k:   代表目前进展到第几层（1——d+2）
a：  第k-1层的节点下标
d：  代表小数位数
T:   每只蚂蚁的路径矩阵
n:   第n只蚂蚁
（注意这个函数竟然和蚂蚁的编号无关，换句话说，
不管哪只蚂蚁他选择路径依据的是同一个函数）
%}
function b=select_k_city(tao,Q0,n,k,T,d)
%test5
%d=7;
%tao=zeros(d+2,10,10);
a=T(n,k-1)+1;%第n只蚂蚁第k-1层的所在城市编码（注意：城市编码从1——10）
q=rand;
%Q0=0.8;
num=100;%轮盘赌次数
P=[];
Select=[];
if q<Q0
    b=find(tao(k,a,:)==max(tao(k,a,:)));  
else
    P=tao(k,a,:)/sum(tao(k,a,:));
    
    Select=Roulette(P,num);

    [row,col,len]=size(tao);
    for i=1:len
     Ps(i)=(sum(Select==i)/num);
    end
    
    b=find(Ps==max(Ps));
end