%{
6) 每只蚂蚁选择城市后都立即按公式(3) 执行局部更新规则;
需要输入rou,tao0,T,tao
%}
%k    第k层
%tao  两层之间的信息素
%rou  比例
%tao0 剩余信息素0.01
%T    每一只蚂蚁下一步到达的城市.第一维是第几只蚂蚁、第二位是目前是第几步
%n    第几只蚂蚁
function tao=update_tao(tao,rou,tao0,T,k,n)
%rou=0.8;
%tao0=0.01;
i=T(n,k)+1;
j=T(n,k-1)+1;
tao(k,j,i)=(1-rou)*tao(k,j,i)+rou*tao0;