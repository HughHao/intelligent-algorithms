function Path_plan()
load map.mat
MM=size(G,1);                      % G ����ͼΪ01�������Ϊ1��ʾ�ϰ���
Tau=ones(MM*MM,MM*MM);             % Tau ��ʼ��Ϣ�ؾ���
Tau=8.*Tau;
K=100;                             %����������ָ���ϳ������ٲ���
M=50;                              %���ϸ���
S=3 ;                              %���·������ʼ��
E=398;                             %���·����Ŀ�ĵ�
Alpha=1;                           % Alpha ������Ϣ����Ҫ�̶ȵĲ���
Beta=7;                            % Beta ��������ʽ������Ҫ�̶ȵĲ���
Rho=0.3 ;                          % Rho ��Ϣ������ϵ��
Q=1;                               % Q ��Ϣ������ǿ��ϵ��
minkl=inf;                         %��ʼ���·�߳���
mink=0;                            %��ʼ��ʹ·����̵���С��������
minl=0;                            %��ʼ��ʹ·����̵���С����
D=G2D(G);
N=size(D,1);                       %N��ʾ����Ĺ�ģ�����ظ�����400
a=1;                               %С�������صı߳�
%�ҳ��յ����ڷ�������ġ�����������������������
Ex=a*(mod(E,MM)-0.5);               %��ֹ�������
if Ex==-0.5
    Ex=MM-0.5;
end
Ey=a*(MM+0.5-ceil(E/MM));            %��ֹ��������ceil ������������Ĵ�����Բ��
Eta=zeros(N);                        %����ʽ��Ϣ������������֮��ת�Ƶ�����

%��������������������������������������������������������������������
%��������ʽ��Ϣ����
for i=1:N                               %N��ʾ����Ĺ�ģ�����ظ�����400
    ix=a*(mod(i,MM)-0.5);               %�õ�ĺ�����
    if ix==-0.5
        ix=MM-0.5;
    end
    iy=a*(MM+0.5-ceil(i/MM));           %ceil ������������Ĵ�����Բ��
    if i~=E                             %E���յ�
        Eta(i)=1/((ix-Ex)^2+(iy-Ey)^2)^0.5; %�õ㵽�յ�ľ���ĵ���
    else
        Eta(i)=100;
    end
end                   %K�ǵ�������100��M����������50
%����������������������������������������������������������������������������
%����ϸ���ṹ�;���洢ÿֻ�����ڸ����е�·���Լ����г��ȡ�������������������
ROUTES=cell(K,M);     %��ϸ���ṹ�洢ÿһ����ÿһֻ���ϵ�����·��
PL=zeros(K,M);         %�þ���洢ÿһ����ÿһֻ���ϵ�����·�߳���
%V0=zeros(K,M);
GASELEC=zeros(K,M);
minGASELEC=inf;
%��������������������������������������������������������������������������������
%����K��������ʳ���ÿ���ɳ�Mֻ����
for k=1:K %K=100
    for m=1:M %M=50
        %״̬��ʼ����ÿֻ���϶��Ǵ�S�㿪ʼ
        W=S;                  %��ǰ�ڵ��ʼ��Ϊ��ʼ��when=start
        Path=S;               %����·�߳�ʼ��
        PLkm=0;               %����·�߳��ȳ�ʼ��pathlongth
        TABUkm=ones(N);       %���ɱ��ʼ��400*400����1
        TABUkm(S)=0;          %�Ѿ��ڳ�ʼ���ˣ����Ҫ�ų�
        DD=D;                 %�ڽӾ����ʼ�������Ⱦ���
        Elec=0; % ��ʼ��������
        Gas=0;  % ��ʼ��������
        v0(1)=10;   % ��ʼ�ٶ�
        t=1;
        v_max=80;% ����ٶ�
        %��һ������ǰ���Ľڵ�
        DW=DD(W,:);   %W��ǰ�㣬DW�Ǿ��뼯�ϣ���W�㰤�ŵĵ������㣬����0
        DW1=find(DW); %���벻��0�Ľڵ�λ�úţ��ڷ��㣨Ҳ���ǰ��ŷ��ϰ���ĵ�ѡ��
        for j=1:length(DW1) %�հ�λ�õ�����
            if TABUkm(DW1(j))==0 %�����λ���Ѿ��ڽ��ɱ���
                DW(j)=0; %��ô�ýڵ㵽��ǰ�ڵ����Ϊ0���Ӷ�ѡȡ�����
            end
        end
        LJD=find(DW);%�ҳ��ڽӿհ׵�ڵ�
        Len_LJD=length(LJD);%��ѡ�ڵ�ĸ���
        
        %����δ����ʳ�������������ͬ������ʳֹͣ
        while W~=E&&Len_LJD>=1 %���統ǰ�ڵ㲻���յ���߿�ѡ�ڵ�������ڵ���1
            %ת�ֶķ�ѡ����һ����ô�ߡ�����������������������������������+++++++++
            PP=zeros(Len_LJD);
            for i=1:Len_LJD
                PP(i)=(Tau(W,LJD(i))^Alpha)*((Eta(LJD(i)))^Beta);
            end
            PP=PP/sum(PP);%�������ʷֲ�
            Pcum=cumsum(PP);
            Select=find(Pcum>=rand);
            to_visit=LJD(Select(1)); %Ҫ���ʽڵ�
            %״̬���ºͼ�¼++++++++++++++++++++++++++++++++
            Path=[Path,to_visit];          %·������
            t=t+1;
            PLkm=PLkm+DD(W,to_visit);     %·����������
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
            %%%%% �ܺ�����
            W=to_visit;                   %�����Ƶ���һ���ڵ�
            for kk=1:N %����DD
                if TABUkm(kk)==0         %�Ѿ����ʵĽڵ��ڽ��ɱ���
                    DD(W,kk)=0;              %W��kk����Ϊ0��Ϊ���´�ѡȡ����㣬�����ظ�
                    DD(kk,W)=0;            %��֮һ���������Ѿ��߹���·
                end
            end
            TABUkm(W)=0;                %�ѷ��ʹ��Ľڵ������ɱ�
            %��һ������ǰ���Ľڵ�
            DW=DD(W,:);   %W��ǰ�㣬DW�Ǿ��뼯��
            DW1=find(DW); %���벻��0�����ţ��Ľڵ�λ�úţ������ڵĵ������0
            for j=1:length(DW1) %���ŵ�λ�õ�����
                if TABUkm(DW1(j))==0 %�����λ���Ѿ��ڽ��ɱ���
                    DW(j)=0; %��ô�ýڵ㵽��ǰ�ڵ����Ϊ0
                end
            end
            LJD=find(DW);%�ҳ������0�Ľڵ�
            Len_LJD=length(LJD);%��ѡ�ڵ�ĸ���
        end
        %����ÿһ��ÿһֻ���ϵ���ʳ·�ߺ�·�߳���
        ROUTES{k,m}=Path; %ϸ���ṹ��·�߼��ϣ���¼ÿ��ÿֻ����·��
        if Path(end)==E %�ߵ��յ�ʱ
            PL(k,m)=PLkm; %Pl��·�����ȣ�k�ǵ���������m�����Ϻ�
            GASELEC(k,m)=Gas+Elec;
            v0(t)=0;
%             if PL(k,m)<minkl %��ǰ����ֻ���ϵ�·�߳���<���·������
%                 mink=k;
%                 minl=m;
%                 minkl=PL(k,m); %��С������������С�����������·������
%             end
            if GASELEC(k,m)<minGASELEC
                mink=k;
                minl=m;
                %minkl=PL(k,m);
                minGASELEC=GASELEC(k,m);
            end
        end
        %���һֻ���ϵ�·��Ѱ��
    end
    
    %������������������������������������������������������������������������������
    %������Ϣ��
    Delta_Tau=zeros(N,N);%��������ʼ��
    for m=1:M
        if PL(k,m)
            ROUT=ROUTES{k,m}; %ϸ���ṹ��·�߼��ϣ���¼ÿ��ÿֻ����·��
            TS=length(ROUT)-1;%����
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
    Tau=(1-Rho).*Tau+Delta_Tau;%��Ϣ�ػӷ�һ���֣�������һ����
    
end %���һ�ε���
%������������������������������������������������������������������������������
%��ͼ
plotif=1;%�Ƿ��ͼ�Ŀ��Ʋ���
if plotif==1 %����������
    %GASELEC=100*ones(K,1); %K��100
    for i=1:K
        PLK=GASELEC(i,:); %PL(i,:)�ǵ�i����������·�����ȣ���50��
        Nonzero= find(PLK); %�ҳ���0����
        PLKPLK=PLK(Nonzero); %ȥ����i��������·�����ȷ�0��ֵ
        minPL(i)=min(PLKPLK); %�õ��ô���̳���
    end
    [H,n]=min(minPL);
    Q=H(1);
    disp(['����ܺ�',num2str(Q)]);b=n(1);
    figure(1)
    plot(minPL(:),'r');
    hold on
    grid on
    title('����ܺı仯����');
    xlabel('��������');
    ylabel('����ܺ�');
    text(0,55,['�����������������ܺ� ',num2str(Q),' �ڵ� ',num2str(b),' ���ﵽ']);
    
    %����������������������������������������������������
    %��������������ͼ
    figure(2)
    axis([0,MM,0,MM]) %�����ڰ׷��񡪡�������������
    for i=1:MM %��ͼ�ϰ����λ��
        for j=1:MM
            if G(i,j)==1 %1���ϰ���
                x1=j-1;y1=MM-i; %�ϰ���������½�
                x2=j;y2=MM-i; %���½�
                x3=j;y3=MM-i+1; %���Ͻ�
                x4=j-1;y4=MM-i+1; %���Ͻ�
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]); %���0.2�Ϻ�
                hold on
            elseif G(i,j)==2%��ͣ
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.5,0.5,0.51]); %���
                hold on
            else
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]); %���
            end
        end
    end
    hold on
    title('�����˶��켣');
    xlabel('����x');
    ylabel('����y');
    plot(2.5,19.5,'go')
    text(2.5,19.5,'���')
    plot(17.5,0.5,'ro')
    text(17.5,0.5,'�յ�')
    %��������·�������·�߹켣������������������������
    ROUT=ROUTES{mink,minl};
    LENROUT=length(ROUT); %���·������
    Rx=ROUT; %������ĺ�����
    Ry=ROUT;
    for ii=1:LENROUT %�ҳ�·�����ڷ��������
        Rx(ii)=a*(mod(ROUT(ii),MM)-0.5);
        if Rx(ii)==-0.5
            Rx(ii)=MM-0.5;
        end
        Ry(ii)=a*(MM+0.5-ceil(ROUT(ii)/MM));
    end
    plot(Rx,Ry,'LineWidth',3)
end
%����������������������������������������������������������2����3���ָ���
%����������������ͼ
plotif2=0;
if plotif2==0
    figure(3)
    %�����ڰ׸񡪡�������������������������������������
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
            elseif G(i,j)==2%��ͣ
                x1=j-1;y1=MM-i;
                x2=j;y2=MM-i;
                x3=j;y3=MM-i+1;
                x4=j-1;y4=MM-i+1;
                fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.5,0.5,0.51]); %���
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
    title('�����˶��켣');
    xlabel('����x');
    ylabel('����y');
    plot(2.5,19.5,'go')
    text(2.5,19.5,'���')
    plot(17.5,0.5,'ro')
    text(17.5,0.5,'�յ�')
    %����·��ͼ��������������������������������
    for k=1:K
        PLK=PL(k,:);
        Nonzero= find(PLK); %�ҳ���0λ��
        PLKPLK=PLK(Nonzero);
        minPLK=min(PLKPLK);
        pos=find(PLK==minPLK); %�ҳ��ô����·���ĸ���λ��
        m=pos(1); %��һ��λ��Ҳ������һk���������ҵ����·�ߵ����Ϻ�m
        ROUT=ROUTES{k,m}; %�ҵ�k��m���ϵ�·��
        LENROUT=length(ROUT); %��������
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
    title('�ٶȱ仯����');
    xlabel('ʱ�䣨s)');
    ylabel('�ٶȣ�km/h)');
%����������������������������������������������������
function D=G2D(G) %������⺯��
l=size(G,1);
D=zeros(l*l,l*l); %ÿ���㵽������ľ����ʼ��0
for i=1:l
    for j=1:l
        if G(i,j)==0 %(i,j)���ϰ���
            for m=1:l
                for n=1:l
                    if G(m,n)==0 %(m,n)���ϰ���(i,j)�ͣ�m,n���ֱ���������������
                        im=abs(i-m);jn=abs(j-n); %����(i,j)��(m,n)������������ֵ
                        %���ҽ�������������ʱ�������
                        if im+jn==1||(im==1&&jn==1) %��������ֵ֮�͵���1����2
                            D((i-1)*l+j,(m-1)*l+n)=(im+jn)^0.5; %���߾���Ϊ1���߸���2
                        end
                    end
                end
            end
        end
        if G(i,j)==2 %ֹͣ�ź�
            for m=1:l
                for n=1:l
                    if G(m,n)==0 %(m,n)���ϰ���(i,j)�ͣ�m,n���ֱ���������������
                        im=abs(i-m);
                        jn=abs(j-n); %����(i,j)��(m,n)������������ֵ
                        %���ҽ�������������ʱ�������
                        if im+jn==1||(im==1&&jn==1) %��������ֵ֮�͵���1����2
                            D((i-1)*l+j,(m-1)*l+n)=2; %���߾���Ϊ1���߸���2
                        end
                    end
                end
            end
        end
        if G(i,j)==0 %ֹͣ�ź�
            for m=1:l
                for n=1:l
                    if G(m,n)==2 %(m,n)���ϰ���(i,j)�ͣ�m,n���ֱ���������������
                        im=abs(i-m);
                        jn=abs(j-n); %����(i,j)��(m,n)������������ֵ
                        %���ҽ�������������ʱ�������
                        if im+jn==1||(im==1&&jn==1) %��������ֵ֮�͵���1����2
                            D((i-1)*l+j,(m-1)*l+n)=2; %���߾���Ϊ1���߸���2
                        end
                    end
                end
            end
        end
    end
end