%% 初始化
%rand 产生每个单元为{0，1}行数为popsize,列数为chromlength的矩阵
%rand 对矩阵的每个单元进行元整，这样产生的初始种群
function pop=fu(popsize,chromlength)
%rand 随机产生每个单元为。。。
%rand对矩阵的每个单元进行元整。。。
pop=round(rand(popsize,chromlength));



