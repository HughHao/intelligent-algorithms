%% 计算目标函数值
%实现计算
function [objvalue]=calobjvalue(pop)
temp1=decodechrom(pop,1,10);%将pop每行转换为十进制
x=temp1*10/1023;%将二值域中数转换为变量域的数
objvalue=10*sin(5*x)+7*cos(4*x);%计算目标函数值