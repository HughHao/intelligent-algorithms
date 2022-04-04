%{
test7
%7) 根据公式(4) ～ (6) 评选出最优蚂蚁并执行全局更新规则;
%共有n只蚂蚁
%tao层与层之间的信息素
d是小数点后精确到d位

%X所有蚂蚁的具体数值数组
%idx最小的蚂蚁编码，也就是数组下标
%alpha是一常量
%}
function [tao,X,idx]=update_all_tao(tao,n,alpha,d,T)
X=zeros(n,1);
for j=1:n
    for i=2:d+1
        X(j)=X(j)+T(j,i)*(10^(1-i));
    end
end

%选择出蚂蚁的最大值
F=[];
for i=1:n
    f=myfunction(X(i));
    F=[F f];
end

cur_max=max(F);
idx=find(F==cur_max);

f=F(idx(1));
%蚂蚁经过路径全局性更新
%alpha=0.8;
for k=2:d+1
    i=T(idx,k-1)+1;
    j=T(idx,k)+1;
    tao(k,i,j)=(1-alpha)*tao(k,i,j)+alpha*(1./f);
end