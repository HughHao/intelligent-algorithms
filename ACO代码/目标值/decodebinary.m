%% 计算目标函数值
%产生[2^n 2^(n-1)...1]的行向量，然后求和，讲二进制转化为十进制
function pop2=decodebinary(pop)
%求pop行和列数
[px,py]=size(pop);
for i = 1:py
    pop1(:,i)=2.^(py-1).*pop(:,i);
end
%求pop1的没行之和
pop2=sum(pop1,2);