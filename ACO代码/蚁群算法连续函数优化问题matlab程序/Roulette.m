%{
���̶ķ�����
%}
function Select=Roulette(P,num)
%:�����̶Ĳ���ѡ����һ��,����num�����̶Ľ��

%:��һ�����̶ķ���,���Ⱥܵ�,
% m = length(P); 
% Select = zeros(1,num);
% for i=1:num
%     Select(i) = m;% ��ʼ��Ϊ���һ��
%     for j=1:m %:������ѡ��
%         if P(j)>rand()
%             Select(i)=j;
%             break;
%         end
%     end
% end

%:�ڶ������̶ķ���,���Ƚϸ�
m = length(P);
Select = zeros(1,num);
r = rand(1,num);
for i=1:num
    sumP = 0;
    j = ceil(m*rand); %����1~m֮����������
    while sumP < r(i)
        sumP = sumP + P(mod(j-1,m)+1);
        j = j+1;
    end
    %Select(i) = mod(j-1,m)+1-1;
    Select(i) = mod(j-2,m)+1;
end

