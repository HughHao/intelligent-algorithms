%% 车间布局遗传禁忌搜索算法仿真主界面
%% 第一步：设置问题实例
Li=[38;16;30;40;48;32;46];%车间长度
Wi=[28;36;16;18;24;28;16];%车间宽度
%每单位距离每单位物流量的物料搬运费用
Pij=[
0,2,3,2,5,4,4;
0,0,5,2,2,3,4;
0,0,0,1,5,4,3;
0,0,0,0,1,5,1;
0,0,0,0,0,4,5;
0,0,0,0,0,0,1;
0,0,0,0,0,0,0
];
%物料搬运的频率
F=[
0,2,2,1,0,2,1;
0,0,2,1,1,2,2;
0,0,0,2,1,2,1;
0,0,0,0,2,1,2;
0,0,0,0,0,2,1;
0,0,0,0,0,0,1;
0,0,0,0,0,0,0
];
%物流量
Q=[
0,10,6,8,4,6,1;
0, 0,3,2,5,4,4;
0, 0,0,6,8,6,5;
0, 0,0,0,5,8,1;
0, 0,0,0,0,8,1;
0, 0,0,0,0,0,5;
0, 0,0,0,0,0,0
];
%物料搬运速率
V=[
0,4,4,4,4,4,4;
0,0,2,2,2,2,2;
0,0,0,2,2,2,2;
0,0,0,0,3,3,3;
0,0,0,0,0,3,3;
0,0,0,0,0,0,2;
0,0,0,0,0,0,0
];
L=200;%矩形区域的长度，x轴
W=120;%矩形区域的宽度，y轴
minDX=10;%各车间的最小水平间距
minDY=10;%各车间的最小垂直间距
minDS=10;%各车间到区域边界的最小距离
%%
pop_size=400;
max_gen=500;
Pm=0.3;
kc=0.5;
kt=0.5;
PLambda=1000;
PK=1000;
n=size(Pij,1);
LB=zeros(2*n,1);
UB=zeros(2*n,1);
for i=1:n
    LB(2*i-1)=0.5*Li(i)+minDS;
    LB(2*i)=0.5*Wi(i)+minDS;
    UB(2*i-1)=L-0.5*Li(i)-minDS;
    UB(2*i)=W-0.5*Wi(i)-minDS;
end
%% 调用遗传算法
figure(3)
[BESTX,BESTY,ALLX,ALLY]=GAUCP2(max_gen,pop_size,Pm,LB,UB,L,W,Li,Wi,Pij,F,Q,V,minDX,minDY,kc,kt,PLambda,PK);
X=BESTX{max_gen};
disp('遗传算法输出的最优结果为');
disp(X);
figure(4)
PlotFigure(X,Li,Wi,L,W);
function [BESTX,BESTY,ALLX,ALLY]=GAUCP2(K,N,Pm,LB,UB,PL,PW,PLi,PWi,PP,PF,PQ,PV,PminDX,PminDY,Pkc,Pkt,PLambda,PK)
%% 此函数实现遗传算法，用于车间布局优化
%% 输入参数列表
% K        迭代次数
% N        种群规模，要求是偶数
% Pm       变异概率
% LB       决策变量的下界，M×1的向量
% UB       决策变量的上界，M×1的向量
%% 输出参数列表
% BESTX    K×1细胞结构，每一个元素是M×1向量，记录每一代的最优个体
% BESTY    K×1矩阵，记录每一代的最优个体的评价函数值
% ALLX     K×1细胞结构，每一个元素是M×N矩阵，记录全部个体
% ALLY     K×N矩阵，记录全部个体的评价函数值
%% 第一步：
M=length(LB);%决策变量的个数LB=zeros(2*n,1)
%种群初始化，每一列是一个样本
farm=zeros(M,N);
for i=1:M
    x=unifrnd(LB(i),UB(i),1,N);
    farm(i,:)=x;
end
%输出变量初始化
ALLX=cell(K,1);%细胞结构，每一个元素是M×N矩阵，记录每一代的个体
ALLY=zeros(K,N);%K×N矩阵，记录每一代评价函数值
BESTX=cell(K,1);%细胞结构，每一个元素是M×1向量，记录每一代的最优个体
BESTY=zeros(K,1);%K×1矩阵，记录每一代的最优个体的评价函数值
k=1;%迭代计数器初始化
%% 第二步：迭代过程
while k<=K
%% 以下是交叉过程
    newfarm=zeros(M,2*N);
    Ser=randperm(N);%两两随机配对的配对表
    A=farm(:,Ser(1));
    B=farm(:,Ser(2));
    P0=unidrnd(M-1);
    a=[A(1:P0,:);B((P0+1):end,:)];%产生子代a
    b=[B(1:P0,:);A((P0+1):end,:)];%产生子代b
    newfarm(:,2*N-1)=a;%加入子代种群
    newfarm(:,2*N)=b;   
    for i=1:(N-1)
        A=farm(:,Ser(i));
        B=farm(:,Ser(i+1));
        P0=unidrnd(M-1);
        a=[A(1:P0,:);B((P0+1):end,:)];
        b=[B(1:P0,:);A((P0+1):end,:)];
        newfarm(:,2*i-1)=a;
        newfarm(:,2*i)=b;
    end   
    FARM=[farm,newfarm];   
%% 选择复制
    SER=randperm(3*N);
    FITNESS=zeros(1,3*N);
    fitness=zeros(1,N);
    for i=1:(3*N)
        Beta=FARM(:,i);
        SE=FIT(Beta,PL,PW,PLi,PWi,PP,PF,PQ,PV,PminDX,PminDY,Pkc,Pkt,PLambda,PK);
        FITNESS(i)=SE;
    end   
    for i=1:N
        f1=FITNESS(SER(3*i-2));
        f2=FITNESS(SER(3*i-1));
        f3=FITNESS(SER(3*i));
        if f1<=f2&&f1<=f3
            farm(:,i)=FARM(:,SER(3*i-2));
            fitness(:,i)=FITNESS(:,SER(3*i-2));
        elseif f2<=f1&&f2<=f3
            farm(:,i)=FARM(:,SER(3*i-1));
            fitness(:,i)=FITNESS(:,SER(3*i-1));
        else
            farm(:,i)=FARM(:,SER(3*i));
            fitness(:,i)=FITNESS(:,SER(3*i));
        end
    end   
    %% 记录最佳个体和收敛曲线
    X=farm;
    Y=fitness;
    ALLX{k}=X;
    ALLY(k,:)=Y;
    minY=min(Y);
    pos=find(Y==minY);
    BESTX{k}=X(:,pos(1));
    BESTY(k)=minY;   
    %% 变异
    for i=1:N
        if Pm>rand&&pos(1)~=i
            AA=farm(:,i);
            BB=GaussMutation(AA,LB,UB);
            for j=1:M
                BB(j,1)=unifrnd(LB(j),UB(j),1,1);
            end
            farm(:,i)=BB;
        end
    end
    disp(k);
    k=k+1;
end
%% 绘图
BESTY2=BESTY;
BESTX2=BESTX;
for k=1:K
    TempY=BESTY(1:k);
    minTempY=min(TempY);
    posY=find(TempY==minTempY);
    BESTY2(k)=minTempY;
    BESTX2{k}=BESTX{posY(1)};
end
BESTY=BESTY2;
BESTX=BESTX2;
MeanBESTY=mean(ALLY');
plot(-BESTY,'-ks','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',1)
hold on
plot(-MeanBESTY,'-ro','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',1)
ylabel('Fitness');
xlabel('Iterations');
legend('Best Fitness','Average Fitness','FontName','Times New Roman','FontSize',10)
grid on
end