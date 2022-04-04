%{
轮盘赌法程序
%}
function Select=Roulette(P,num)
%:按轮盘赌策略选择下一点,返回num次轮盘赌结果

%:第一种轮盘赌方法,精度很低,
% m = length(P); 
% Select = zeros(1,num);
% for i=1:num
%     Select(i) = m;% 初始化为最后一个
%     for j=1:m %:按概率选择
%         if P(j)>rand()
%             Select(i)=j;
%             break;
%         end
%     end
% end

%:第二种轮盘赌方法,精度较高
m = length(P);
Select = zeros(1,num);
r = rand(1,num);
for i=1:num
    sumP = 0;
    j = ceil(m*rand); %产生1~m之间的随机整数
    while sumP < r(i)
        sumP = sumP + P(mod(j-1,m)+1);
        j = j+1;
    end
    %Select(i) = mod(j-1,m)+1-1;
    Select(i) = mod(j-2,m)+1;
end

