%{
test7
%7) ���ݹ�ʽ(4) �� (6) ��ѡ���������ϲ�ִ��ȫ�ָ��¹���;
%����nֻ����
%tao�����֮�����Ϣ��
d��С�����ȷ��dλ

%X�������ϵľ�����ֵ����
%idx��С�����ϱ��룬Ҳ���������±�
%alpha��һ����
%}
function [tao,X,idx]=update_all_tao(tao,n,alpha,d,T)
X=zeros(n,1);
for j=1:n
    for i=2:d+1
        X(j)=X(j)+T(j,i)*(10^(1-i));
    end
end

%ѡ������ϵ����ֵ
F=[];
for i=1:n
    f=myfunction(X(i));
    F=[F f];
end

cur_max=max(F);
idx=find(F==cur_max);

f=F(idx(1));
%���Ͼ���·��ȫ���Ը���
%alpha=0.8;
for k=2:d+1
    i=T(idx,k-1)+1;
    j=T(idx,k)+1;
    tao(k,i,j)=(1-alpha)*tao(k,i,j)+alpha*(1./f);
end