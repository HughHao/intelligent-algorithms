%test_lianxu_function
%��һ������һ��ֻ��0
%�м�㶼��0-9  len
%T���ϵ�·��������ʼ���ǵ�һ���0���յ��ǵ�d+2���0
%taoȫ��d+2�����Ϣ�ؾ��󣬴���ÿһ����ÿһ�ڵ�ľ�������
len=10;
tao0=0.01;
d=7;
tao=ones(d+2,len,len);
tao=tao*tao0;%�������Ϲ���
rou=0.8;
N0=20;
T=zeros(N0,d+2);%N0ֻ���ϣ�d+2�㣬�����ʼֵΪ0��Ҳ����˵�������ϵ�һ��εĶ��ӵ�һ�����г���
Q0=0.8;
d=7;
iteration=100;
N0=20;
alpha=0.8;

max_result=0;
max_route=[];
max_x=0;
for i=1:iteration%��������������1000�Σ�Ҳ�������������߹�1000��
    %a=1;%��һ��������±�1��ʼ
    for j=2:d+1%�������ӳ�ʼ���һ��0�����������һ���0����
        %b=1;
        for k=1:N0%���ϱ���
            b=select_k_city(tao,Q0,k,j,T,d);
            %b=select_k_city(tao,Q0,j,a,d);%�����kֻ���ϣ���j��ڵ�λ��
            T(k,j)=b(1)-1;%T���鱣����ǵ�k���ϵ�j��ľ�����ֵ��Ӧ��Ϊ���б�ż�1
            
            %ִ�оֲ����³���֮�����Ϣ��
            %k    ��k��
            %tao  ����֮�����Ϣ��
            %rou  ����
            %tao0 ʣ����Ϣ��0.01
            %T    ÿһֻ������һ������ĳ���.��һά�ǵڼ�ֻ���ϡ��ڶ�λ��Ŀǰ�ǵڼ���
            %n    �ڼ�ֻ����
            tao=update_tao(tao,rou,tao0,T,j,k);
        end
        %a=b;
    end
    
    %����������ȫ�ֶ�����һ��󣬽���ȫ���Ը���
    %7) ���ݹ�ʽ(4) �� (6) ��ѡ���������ϲ�ִ��ȫ�ָ��¹���;
    %����nֻ����
    %tao�����֮�����Ϣ��
    %d��С�����ȷ��dλ

    %X�������ϵľ�����ֵ����
    %idx��С�����ϱ��룬Ҳ���������±�
    %alpha��һ����
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
