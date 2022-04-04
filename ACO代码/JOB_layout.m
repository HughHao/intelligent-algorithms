%% ���䲼���Ŵ����������㷨����������
%% ��һ������������ʵ��
Li=[38;16;30;40;48;32;46];%���䳤��
Wi=[28;36;16;18;24;28;16];%������
%ÿ��λ����ÿ��λ�����������ϰ��˷���
Pij=[
0,2,3,2,5,4,4;
0,0,5,2,2,3,4;
0,0,0,1,5,4,3;
0,0,0,0,1,5,1;
0,0,0,0,0,4,5;
0,0,0,0,0,0,1;
0,0,0,0,0,0,0
];
%���ϰ��˵�Ƶ��
F=[
0,2,2,1,0,2,1;
0,0,2,1,1,2,2;
0,0,0,2,1,2,1;
0,0,0,0,2,1,2;
0,0,0,0,0,2,1;
0,0,0,0,0,0,1;
0,0,0,0,0,0,0
];
%������
Q=[
0,10,6,8,4,6,1;
0, 0,3,2,5,4,4;
0, 0,0,6,8,6,5;
0, 0,0,0,5,8,1;
0, 0,0,0,0,8,1;
0, 0,0,0,0,0,5;
0, 0,0,0,0,0,0
];
%���ϰ�������
V=[
0,4,4,4,4,4,4;
0,0,2,2,2,2,2;
0,0,0,2,2,2,2;
0,0,0,0,3,3,3;
0,0,0,0,0,3,3;
0,0,0,0,0,0,2;
0,0,0,0,0,0,0
];
L=200;%��������ĳ��ȣ�x��
W=120;%��������Ŀ�ȣ�y��
minDX=10;%���������Сˮƽ���
minDY=10;%���������С��ֱ���
minDS=10;%�����䵽����߽����С����
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
%% �����Ŵ��㷨
figure(3)
[BESTX,BESTY,ALLX,ALLY]=GAUCP2(max_gen,pop_size,Pm,LB,UB,L,W,Li,Wi,Pij,F,Q,V,minDX,minDY,kc,kt,PLambda,PK);
X=BESTX{max_gen};
disp('�Ŵ��㷨��������Ž��Ϊ');
disp(X);
figure(4)
PlotFigure(X,Li,Wi,L,W);
function [BESTX,BESTY,ALLX,ALLY]=GAUCP2(K,N,Pm,LB,UB,PL,PW,PLi,PWi,PP,PF,PQ,PV,PminDX,PminDY,Pkc,Pkt,PLambda,PK)
%% �˺���ʵ���Ŵ��㷨�����ڳ��䲼���Ż�
%% ��������б�
% K        ��������
% N        ��Ⱥ��ģ��Ҫ����ż��
% Pm       �������
% LB       ���߱������½磬M��1������
% UB       ���߱������Ͻ磬M��1������
%% ��������б�
% BESTX    K��1ϸ���ṹ��ÿһ��Ԫ����M��1��������¼ÿһ�������Ÿ���
% BESTY    K��1���󣬼�¼ÿһ�������Ÿ�������ۺ���ֵ
% ALLX     K��1ϸ���ṹ��ÿһ��Ԫ����M��N���󣬼�¼ȫ������
% ALLY     K��N���󣬼�¼ȫ����������ۺ���ֵ
%% ��һ����
M=length(LB);%���߱����ĸ���LB=zeros(2*n,1)
%��Ⱥ��ʼ����ÿһ����һ������
farm=zeros(M,N);
for i=1:M
    x=unifrnd(LB(i),UB(i),1,N);
    farm(i,:)=x;
end
%���������ʼ��
ALLX=cell(K,1);%ϸ���ṹ��ÿһ��Ԫ����M��N���󣬼�¼ÿһ���ĸ���
ALLY=zeros(K,N);%K��N���󣬼�¼ÿһ�����ۺ���ֵ
BESTX=cell(K,1);%ϸ���ṹ��ÿһ��Ԫ����M��1��������¼ÿһ�������Ÿ���
BESTY=zeros(K,1);%K��1���󣬼�¼ÿһ�������Ÿ�������ۺ���ֵ
k=1;%������������ʼ��
%% �ڶ�������������
while k<=K
%% �����ǽ������
    newfarm=zeros(M,2*N);
    Ser=randperm(N);%���������Ե���Ա�
    A=farm(:,Ser(1));
    B=farm(:,Ser(2));
    P0=unidrnd(M-1);
    a=[A(1:P0,:);B((P0+1):end,:)];%�����Ӵ�a
    b=[B(1:P0,:);A((P0+1):end,:)];%�����Ӵ�b
    newfarm(:,2*N-1)=a;%�����Ӵ���Ⱥ
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
%% ѡ����
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
    %% ��¼��Ѹ������������
    X=farm;
    Y=fitness;
    ALLX{k}=X;
    ALLY(k,:)=Y;
    minY=min(Y);
    pos=find(Y==minY);
    BESTX{k}=X(:,pos(1));
    BESTY(k)=minY;   
    %% ����
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
%% ��ͼ
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