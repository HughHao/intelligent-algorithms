%% ����Ŀ�꺯��ֵ
%����[2^n 2^(n-1)...1]����������Ȼ����ͣ���������ת��Ϊʮ����
function pop2=decodebinary(pop)
%��pop�к�����
[px,py]=size(pop);
for i = 1:py
    pop1(:,i)=2.^(py-1).*pop(:,i);
end
%��pop1��û��֮��
pop2=sum(pop1,2);