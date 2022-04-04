function Path_plan()
load map.mat
MM=size(G,1);                      % G 地形图为01矩阵，如果为1表示障碍物
Tau=ones(MM*MM,MM*MM);             % Tau 初始信息素矩阵
Tau=8.*Tau;
K=100;                             %迭代次数（指蚂蚁出动多少波）
M=50;                              %蚂蚁个数
S=3 ;                              %最短路径的起始点
E=398;                             %最短路径的目的点
Alpha=1;                           % Alpha 表征信息素重要程度的参数
Beta=7;                            % Beta 表征启发式因子重要程度的参数
Rho=0.3 ;                          % Rho 信息素蒸发系数
Q=1;                               % Q 信息素增加强度系数
minkl=inf;                         %初始最短路线长度
mink=0;                            %初始化使路径最短的最小迭代次数
minl=0;                            %初始化使路径最短的最小蚂蚁
D=G2D(G);
N=size(D,1);                       %N表示问题的规模（象素个数）400
a=1;                               %小方格象素的边长
%找出终点所在方格的中心————————————
Ex=a*(mod(E,MM)-0.5);               %终止点横坐标
if Ex==-0.5
    Ex=MM-0.5;
end
Ey=a*(MM+0.5-ceil(E/MM));            %终止点纵坐标ceil 是向离它最近的大整数圆整
Eta=zeros(N);                        %启发式信息，人与两个点之间转移的期望

%——————————————————————————————————
%以下启发式信息矩阵
for i=1:N                               %N表示问题的规模（象素个数）400
    ix=a*(mod(i,MM)-0.5);               %该点的横坐标
    if ix==-0.5
        ix=MM-0.5;
    end
    iy=a*(MM+0.5-ceil(i/MM));           %ceil 是向离它最近的大整数圆整
    if i~=E                             %E是终点
        Eta(i)=1/((ix-Ex)^2+(iy-Ey)^2)^0.5; %该点到终点的距离的倒数
    else
        Eta(i)=100;
    end
end                   %K是迭代次数100，M是蚂蚁数量50
%——————————————————————————————————————
%定义细胞结构和矩阵存储每只蚂蚁在各代中的路线以及爬行长度——————————
ROUTES=cell(K,M);     %用细胞结构存储每一代的每一只蚂蚁的爬行路线
PL=zeros(K,M);         %用矩阵存储每一代的每一只蚂蚁的爬行路线长度
%V0=zeros(K,M);
GASELEC=zeros(K,M);
minGASELEC=inf;
%————————————————————————————————————————
%启动K轮蚂蚁觅食活动，每轮派出M只蚂蚁
for k=1:K %K=100
    for m=1:M %M=50
        %状态初始化，每只蚂蚁都是从S点开始
        W=S;                  %当前节点初始化为起始点when=start
        Path=S;               %爬行路线初始化
        PLkm=0;               %爬行路线长度初始化pathlongth
        TABUkm=ones(N);       %禁忌表初始化400*400矩阵1
        TABUkm(S)=0;          %已经在初始点了，因此要排除
        DD=D;                 %邻接矩阵初始化，长度矩阵
        Elec=0; % 初始电量消耗
        Gas=0;  % 初始油量消耗
        v0(1)=10;   % 起始速度
        t=1;
        v_max=80;% 最大速度
        %下一步可以前往的节点
        DW=DD(W,:);   %W当前点，DW是距离集合，与W点挨着的点距离非零，其他0
        DW1=find(DW); %距离不是0的节点位置号，在非零（也就是挨着非障碍物）的点选择
        for j=1:length(DW1) %空白位置的数量
            if TABUkm(DW1(j))==0 %如果该位置已经在禁忌表中
                DW(j)=0; %那么该节点到当前节点距离为0，从而选取非零点
            end
        end
        LJD=find(DW);%找出邻接空白点节点
        Len_LJD=length(LJD);%可选节点的个数
        
        %蚂蚁未遇到食物或者陷入死胡同或者觅食停止
        while W~=E&&Len_LJD>=1 %假如当前节点不在终点或者可选节点个数大于等于1
            %转轮赌法选择下一步怎么走——————————————————+++++++++
            PP=zeros(Len_LJD);
            for i=1:Len_LJD
                PP(i)=(Tau(W,LJD(i))^Alpha)*((Eta(LJD(i)))^Beta);
            end
            PP=PP/sum(PP);%建立概率分布
            Pcum=cumsum(PP);
            Select=find(Pcum>=rand);
            to_visit=LJD(Select(1)); %要访问节点
            %状态更新和记录++++++++++++++++++++++++++++++++
            Path=[Path,to_visit];          %路径增加
            t=t+1;
            PLkm=PLkm+DD(W,to_visit);     %路径长度增加
            to_x=a*(mod(to_visit,MM)-0.5);    
            if to_x==-0.5
                to_x=MM-0.5;
            end
            to_y=a*(MM+0.5-ceil(to_visit/MM));
            
            W_x=a*(mod(W,MM)-0.5);    
            if W_x==-0.5
                W_x=MM-0.5;
            end
            W_y=a*(MM+0.5-ceil(W/MM));
            
            if t==2
                v0(t)=v0(t-1)+12*DD(W,to_visit);
                Elec=1;
                Gas=1;
            elseif t>2
                before=Path(t-2);
                before_x=a*(mod(before,MM)-0.5);
                if before_x==-0.5
                    before_x=MM-0.5;
                end
                before_y=a*(MM+0.5-ceil(before/MM));
                if (to_y-W_y)/(to_x-W_x)==(W_y-before_y)/(W_x-before_x)
                    if v0(t-1)+10*DD(W,to_visit)<v_max
                        if t>100
                            v0(t)=DD(W,to_visit);
                        else
                            v0(t)=v0(t-1)+DD(W,to_visit);
                        end
                    else
                        v0(t)=v_max;
                    end
                    Elec=Elec+1;
                    Gas=Gas+1;
                else
                    if v0(t-1)/1.1<20
                        if t>100
                            v0(t)=DD(W,to_visit);
                        else
                            v0(t)=10*60/(t+10)+5;
                        end
                    else
                        v0(t)=v0(t-1)/1.1;
                    end
                    Elec=Elec+0.5;
                    Gas=Gas+2^0.5;
                end
                if DD(Path(t-1),Path(t))==2
                    if t>100
                        v0(t)=DD(W,to_visit);
                    else
                        v0(t)=10*60/(t+10)+1;
                    end
                end
            end
            %%%%% 能耗增加
            W=to_visit;                   %蚂蚁移到下一个节点
            for kk=1:N %更新DD
                if TABUkm(kk)==0         %已经访问的节点在禁忌表中
                    DD(W,kk)=0;              %W到kk距离为0，为了下次选取非零点，不能重复
                    DD(kk,W)=0;            %反之一样，不走已经走过的路
                end
            end
            TABUkm(W)=0;                %已访问过的节点加入禁忌表
            %下一步可以前往的节点
            DW=DD(W,:);   %W当前点，DW是距离集合
            DW1=find(DW); %距离不是0（挨着）的节点位置号，不相邻的点距离设0
            for j=1:length(DW1) %挨着点位置的数量
                if TABUkm(DW1(j))==0 %如果该位置已经在禁忌表中
                    DW(j)=0; %那么该节点到当前节点距离为0
                end
            end
            LJD=find(DW);%找出距离非0的节点
            Len_LJD=length(LJD);%可选节点的个数
        end
        %记下每一代每一只蚂蚁的觅食路线和路线长度
        ROUTES{k,m}=Path; %细胞结构，路线集合，记录每代每只蚂蚁路线
        if Path(end)==E %走到终点时
            PL(k,m)=PLkm; %Pl是路径长度，k是迭代次数，m是蚂蚁号
            GASELEC(k,m)=Gas+Elec;
            v0(t)=0;
%             if PL(k,m)<minkl %当前代当只蚂蚁的路线长度<最短路径长度
%                 mink=k;
%                 minl=m;
%                 minkl=PL(k,m); %最小迭代次数，最小蚂蚁数，最短路径更新
%             end
            if GASELEC(k,m)<minGASELEC
                mink=k;
                minl=m;
                %minkl=PL(k,m);
                minGASELEC=GASELEC(k,m);
            end
        end
        %完成一只蚂蚁的路线寻优
    end
    
    %———————————————————————————————————————
    %更新信息素
    Delta_Tau=zeros(N,N);%更新量初始化
    for m=1:M
        if PL(k,m)
            ROUT=ROUTES{k,m}; %细胞结构，路线集合，记录每代每只蚂蚁路线
            TS=length(ROUT)-1;%跳数
            %PL_km=PL(k,m);
            GAS_ELEC=GASELEC(k,m);
            for s=1:TS
                x=ROUT(s);
                y=ROUT(s+1);
                Delta_Tau(x,y)=Delta_Tau(x,y)+Q/GAS_ELEC;
                Delta_Tau(y,x)=Delta_Tau(y,x)+Q/GAS_ELEC;
            end
        end
    end
    Tau=(1-Rho).*Tau+Delta_Tau;%信息素挥发一部分，新增加一部分
    
end %完成一次迭代
%———————————————————————————————————————
%绘图
plotif=1;%是否绘图的控制参数
if plotif==1 %绘收敛曲线
    %GASELEC=100*ones(K,1); %K是100
    for i=1:K
        PLK=GASELEC(i,:); %PL(i,:)是第i代所有蚂蚁路径长度，共50个
        Nonzero= find(PLK); %找出非0蚂蚁
        PLKPLK=PLK(Nonzero); %去除第i代蚂蚁中路径长度非0的值
        minPL(i)=min(PLKPLK); %得到该代最短长度
    end
    [H,n]=min(minPL);
    Q=H(1);
    disp(['最低能耗',num2str(Q)]);b=n(1);
    figure(1)
    plot(minPL(:),'r');
    hold on
    grid on
    title('最低能耗变化趋势');
    xlabel('迭代次数');
    ylabel('最低能耗');
    text(0,55,['最终搜索结果：最低能耗 ',num2str(Q),' 在第 ',num2str(b),' 代达到']);
    
    %——————————————————————————
    %绘最优蚂蚁爬行图
    figure(2)
    axis([0,MM,0,MM]) %画出黑白方格————————
    for i=1:MM %绘图障碍物的位置
        for j=1:MM
            if G(i,j)==1 %1是障碍物
                x1=j-1;y1=MM-i; %障碍物起点左下角
                x2=j;y2=MM-i; %右下角
                x3=j;y3=MM-i+1; %右上角
                x4=j-1;y4=MM-i+1; %左上角
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]); %填充0.2较黑
                hold on
            elseif G(i,j)==2%暂停
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.5,0.5,0.51]); %填充
                hold on
            else
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]); %填充
            end
        end
    end
    hold on
    title('汽车运动轨迹');
    xlabel('坐标x');
    ylabel('坐标y');
    plot(2.5,19.5,'go')
    text(2.5,19.5,'起点')
    plot(17.5,0.5,'ro')
    text(17.5,0.5,'终点')
    %画出所有路线中最短路线轨迹————————————
    ROUT=ROUTES{mink,minl};
    LENROUT=length(ROUT); %最短路径长度
    Rx=ROUT; %经过点的横坐标
    Ry=ROUT;
    for ii=1:LENROUT %找出路径所在方格的中心
        Rx(ii)=a*(mod(ROUT(ii),MM)-0.5);
        if Rx(ii)==-0.5
            Rx(ii)=MM-0.5;
        end
        Ry(ii)=a*(MM+0.5-ceil(ROUT(ii)/MM));
    end
    plot(Rx,Ry,'LineWidth',3)
end
%————————————————————————————（2）（3）分割线
%绘各代最短蚂蚁爬行图
plotif2=0;
if plotif2==0
    figure(3)
    %画出黑白格————————————————————
    axis([0,MM,0,MM])
    for i=1:MM
        for j=1:MM
            if G(i,j)==1
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]);
                hold on
            elseif G(i,j)==2%暂停
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.5,0.5,0.51]); %填充
                hold on
            else
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]);
                hold on
            end
        end
    end
    title('汽车运动轨迹');
    xlabel('坐标x');
    ylabel('坐标y');
    plot(2.5,19.5,'go')
    text(2.5,19.5,'起点')
    plot(17.5,0.5,'ro')
    text(17.5,0.5,'终点')
    %画出路线图————————————————
    for k=1:K
        PLK=PL(k,:);
        Nonzero= find(PLK); %找出非0位置
        PLKPLK=PLK(Nonzero);
        minPLK=min(PLKPLK);
        pos=find(PLK==minPLK); %找出该代最短路径的各个位置
        m=pos(1); %第一个位置也即是这一k代中最早找到最短路线的蚂蚁号m
        ROUT=ROUTES{k,m}; %找到k代m蚂蚁的路线
        LENROUT=length(ROUT); %方格数量
        Rx=ROUT;
        Ry=ROUT;
        for ii=1:LENROUT
            Rx(ii)=a*(mod(ROUT(ii),MM)-0.5);
            if Rx(ii)==-0.5
                Rx(ii)=MM-0.5;
            end
            Ry(ii)=a*(MM+0.5-ceil(ROUT(ii)/MM));
        end
        plot(Rx,Ry,'LineWidth',1)
        hold on
    end
end
T=180000;
p=1;
while p<T
    q=fix(p/120)+1;
    if q<=t
        V0(p)=v0(q)+rand*v0(q);
        if V0(p)>v_max
            V0(p)=v_max;
        end
        if v0(q)==0 && q<100
            V0(p-10:p-1)=abs(V0(p-10:p-1)/2-p/20);
        end
        
    end
    p=p+1;
end
lenV=length(V0);
for i=lenV-100:lenV
    V0(i)=V0(i-1)*(lenV-i)/100;
end
V0(1)=0;
V0(2)=1;

for j=3:80
    if rand>0.5
        V0(j)=V0(j-1)-0.2*rand;
    else
        V0(j)=V0(j-1)+2*rand;
    end
end
for j=2:(fix(lenV/80)-1)
    for i=80*(j-1)+1:80*j
        if rand>0.5
            V0(i)=V0(i-1)+2*rand;
        else
            V0(i)=abs(V0(i-1)-2*rand);
        end
    end
end
for i=(fix(lenV/80)-1)*80+1:lenV-50
    if V0(i-1)>20
        if rand>0.5
            V0(i)=V0(i-1)-2*rand;
        else
            V0(i)=V0(i-1)-rand;
        end
    else
        V0(i)=V0(i-1)+rand;
    end
end
for i=lenV-49:lenV-1
    if V0(i)>5
        V0(i)=V0(i-1)-0.5*rand;
    else
        V0(i)=V0(i)+0.01*rand;
    end
end
V0(end)=0;
for i=1:lenV
    if V0(i)>v_max
        V0(i)=v_max+rand;
    end
end
xlswrite('v_t3.xlsx',V0')
t
figure(4)
    plot(V0(:));
    title('速度变化曲线');
    xlabel('时间（s)');
    ylabel('速度（km/h)');
%——————————————————————————
function D=G2D(G) %距离求解函数
l=size(G,1);
D=zeros(l*l,l*l); %每个点到其他点的距离初始化0
for i=1:l
    for j=1:l
        if G(i,j)==0 %(i,j)非障碍物
            for m=1:l
                for n=1:l
                    if G(m,n)==0 %(m,n)非障碍物(i,j)和（m,n）分别代表两个点的坐标
                        im=abs(i-m);jn=abs(j-n); %计算(i,j)与(m,n)横纵坐标差绝对值
                        %当且仅当两个方格挨着时计算距离
                        if im+jn==1||(im==1&&jn==1) %坐标差绝对值之和等于1或者2
                            D((i-1)*l+j,(m-1)*l+n)=(im+jn)^0.5; %两者距离为1或者根号2
                        end
                    end
                end
            end
        end
        if G(i,j)==2 %停止信号
            for m=1:l
                for n=1:l
                    if G(m,n)==0 %(m,n)非障碍物(i,j)和（m,n）分别代表两个点的坐标
                        im=abs(i-m);
                        jn=abs(j-n); %计算(i,j)与(m,n)横纵坐标差绝对值
                        %当且仅当两个方格挨着时计算距离
                        if im+jn==1||(im==1&&jn==1) %坐标差绝对值之和等于1或者2
                            D((i-1)*l+j,(m-1)*l+n)=2; %两者距离为1或者根号2
                        end
                    end
                end
            end
        end
        if G(i,j)==0 %停止信号
            for m=1:l
                for n=1:l
                    if G(m,n)==2 %(m,n)非障碍物(i,j)和（m,n）分别代表两个点的坐标
                        im=abs(i-m);
                        jn=abs(j-n); %计算(i,j)与(m,n)横纵坐标差绝对值
                        %当且仅当两个方格挨着时计算距离
                        if im+jn==1||(im==1&&jn==1) %坐标差绝对值之和等于1或者2
                            D((i-1)*l+j,(m-1)*l+n)=2; %两者距离为1或者根号2
                        end
                    end
                end
            end
        end
    end
end