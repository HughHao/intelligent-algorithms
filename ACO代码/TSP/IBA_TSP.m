function IBA_TSP
 
CityNum=31; %you can choose 10, 30, 50, 75
[dislist,Clist]=tsp(CityNum);%�����嵥�������е�
 
inn=31; %��ʼ��Ⱥ��С
gnmax=500; %������
pc=0.8; %�������
pm=0.8; %�������
 
%������ʼ��Ⱥ
zhongqun=zeros(inn,CityNum);
for i=1:inn
    zhongqun(i,:)=randperm(CityNum);%����������С�������������
end
[~,p]=objf(zhongqun,dislist);%Ŀ����Ӧ��ֵ
 
gn=1;
ymean=zeros(gn,1);
ymax=zeros(gn,1);
xmax=zeros(inn,CityNum);
scnew=zeros(inn,CityNum);
smnew=zeros(inn,CityNum);
while gn<gnmax+1
   for j=1:2:inn
      seln=sel(p);  %ѡ�����
      scro=cro(zhongqun,seln,pc);  %�������
      scnew(j,:)=scro(1,:);
      scnew(j+1,:)=scro(2,:);
      smnew(j,:)=mut(scnew(j,:),pm);  %�������
      smnew(j+1,:)=mut(scnew(j+1,:),pm);
   end
   zhongqun=smnew;  %�������µ���Ⱥ
   [f,p]=objf(zhongqun,dislist);  %��������Ⱥ����Ӧ��
   %��¼��ǰ����ú�ƽ������Ӧ��
   [fmax,nmax]=max(f);
   ymean(gn)=1000/mean(f);
   ymax(gn)=1000/fmax;
   %��¼��ǰ������Ѹ���
   x=zhongqun(nmax,:);
   xmax(gn,:)=x;
   drawTSP(Clist,x,ymax(gn),gn,0);
   gn=gn+1;
end
[min_ymax,jn]=min(ymax);
drawTSP(Clist,xmax(jn,:),min_ymax,jn,1);
 
figure(2);
plot(ymax,'r'); hold on;
plot(ymean,'b');grid on;
title('��������');
legend('���Ž�','ƽ����');
fprintf('�Ŵ��㷨�õ�����̾���:%.2f\n',min_ymax);
fprintf('�Ŵ��㷨�õ������·��');
disp(xmax(jn,:));
end
 
%------------------------------------------------
%����������Ⱥ����Ӧ��
function [f,p]=objf(s,dislist)
 
inn=size(s,1);  %��ȡ��Ⱥ��С
f=zeros(inn,1);
for i=1:inn
   f(i)=CalDist(dislist,s(i,:));  %���㺯��ֵ������Ӧ��
end
f=1000./f'; %ȡ���뵹��
%���ݸ������Ӧ�ȼ����䱻ѡ��ĸ���
fsum=0;
for i=1:inn
   fsum=fsum+f(i)^15;% ����Ӧ��Խ�õĸ��屻ѡ�����Խ��
end
ps=zeros(inn,1);
for i=1:inn
   ps(i)=f(i)^15/fsum;
end
 
%�����ۻ�����
p=zeros(inn,1);
p(1)=ps(1);
for i=2:inn
   p(i)=p(i-1)+ps(i);
end
p=p';
end
 
%--------------------------------------------------
%���ݱ�������ж��Ƿ����
function pcc=pro(pc)
test(1:100)=0;
l=round(100*pc);
test(1:l)=1;
n=round(rand*99)+1;
pcc=test(n);   
end
 
%--------------------------------------------------
%��ѡ�񡱲���
function seln=sel(p)
 
seln=zeros(2,1);
%����Ⱥ��ѡ���������壬��ò�Ҫ����ѡ��ͬһ������
for i=1:2
   r=rand;  %����һ�������
   prand=p-r;
   j=1;
   while prand(j)<0
       j=j+1;
   end
   seln(i)=j; %ѡ�и�������
   if i==2&&j==seln(i-1)    %%����ͬ����ѡһ��
       r=rand;  %����һ�������
       prand=p-r;
       j=1;
       while prand(j)<0
           j=j+1;
       end
   end
end
end
 
%------------------------------------------------
%�����桱����
function scro=cro(s,seln,pc)
 
bn=size(s,2);
pcc=pro(pc);  %���ݽ�����ʾ����Ƿ���н��������1���ǣ�0���
scro(1,:)=s(seln(1),:);
scro(2,:)=s(seln(2),:);
if pcc==1
   c1=round(rand*(bn-2))+1;  %��[1,bn-1]��Χ���������һ������λ
   c2=round(rand*(bn-2))+1;
   chb1=min(c1,c2);
   chb2=max(c1,c2);
   middle=scro(1,chb1+1:chb2);
   scro(1,chb1+1:chb2)=scro(2,chb1+1:chb2);
   scro(2,chb1+1:chb2)=middle;
   for i=1:chb1 %�ƺ�������
       while find(scro(1,chb1+1:chb2)==scro(1,i))
           zhi=find(scro(1,chb1+1:chb2)==scro(1,i));
           y=scro(2,chb1+zhi);
           scro(1,i)=y;
       end
       while find(scro(2,chb1+1:chb2)==scro(2,i))
           zhi=find(scro(2,chb1+1:chb2)==scro(2,i));
           y=scro(1,chb1+zhi);
           scro(2,i)=y;
       end
   end
   for i=chb2+1:bn
       while find(scro(1,1:chb2)==scro(1,i))
           zhi=logical(scro(1,1:chb2)==scro(1,i));
           y=scro(2,zhi);
           scro(1,i)=y;
       end
       while find(scro(2,1:chb2)==scro(2,i))
           zhi=logical(scro(2,1:chb2)==scro(2,i));
           y=scro(1,zhi);
           scro(2,i)=y;
       end
   end
end
end
 
%--------------------------------------------------
%�����족����
function snnew=mut(snew,pm)
 
bn=size(snew,2);
snnew=snew;
 
pmm=pro(pm);  %���ݱ�����ʾ����Ƿ���б��������1���ǣ�0���
if pmm==1
   c1=round(rand*(bn-2))+1;  %��[1,bn-1]��Χ���������һ������λ
   c2=round(rand*(bn-2))+1;
   chb1=min(c1,c2);
   chb2=max(c1,c2);
   x=snew(chb1+1:chb2);
   snnew(chb1+1:chb2)=fliplr(x);
end
end
 
%------------------------------------------------
%����λ������
function [dislist,cityn]=tsp(n)
dislist=zeros(n,n);
if n==10
    city10=[0.4 0.4439;0.2439 0.1463;0.1707 0.2293;0.2293 0.761;0.5171 0.9414;
        0.8732 0.6536;0.6878 0.5219;0.8488 0.3609;0.6683 0.2536;0.6195 0.2634];%10 cities d'=2.691
    for i=1:10
        for j=1:10
            dislist(i,j)=((city10(i,1)-city10(j,1))^2+(city10(i,2)-city10(j,2))^2)^0.5;
        end
    end
    cityn=city10;
end
if n==31
    city30=[1304.0	2312.0;3639.0	1315.0;4177.0	2244.0;3712.0	1399.0;3488.0	1535.0;3326.0	1556.0;3238.0	1229.0;
4196.0	1004.0;4312.0	790.0;4386.0	570.0;3007.0	1970.0;2562.0	1756.0;2788.0	1491.0;2381.0	1676.0;1332.0	695.0;
3715.0	1678.0;3918.0	2179.0;4061.0	2370.0;3780.0	2212.0;3676.0	2578.0;4029.0	2838.0;4263.0	2931.0;
3429.0	1908.0
3507.0	2367.0
3394.0	2643.0
3439.0	3201.0
2935.0	3240.0
3140.0	3550.0
2545.0	2357.0
2778.0	2826.0
2370.0	2975.0
];%30 cities d'=423.741 by D B Fogel
    for i=1:31
        for j=1:31
            dislist(i,j)=((city30(i,1)-city30(j,1))^2+(city30(i,2)-city30(j,2))^2)^0.5;
        end
    end
    cityn=city30;
end
 
if n==50
    city50=[31 32;32 39;40 30;37 69;27 68;37 52;38 46;31 62;30 48;21 47;25 55;16 57;
        17 63;42 41;17 33;25 32;5 64;8 52;12 42;7 38;5 25; 10 77;45 35;42 57;32 22;
        27 23;56 37;52 41;49 49;58 48;57 58;39 10;46 10;59 15;51 21;48 28;52 33;
        58 27;61 33;62 63;20 26;5 6;13 13;21 10;30 15;36 16;62 42;63 69;52 64;43 67];%50 cities d'=427.855 by D B Fogel
    for i=1:50
        for j=1:50
            dislist(i,j)=((city50(i,1)-city50(j,1))^2+(city50(i,2)-city50(j,2))^2)^0.5;
        end
    end
    cityn=city50;
end
 
if n==75
    city75=[48 21;52 26;55 50;50 50;41 46;51 42;55 45;38 33;33 34;45 35;40 37;50 30;
        55 34;54 38;26 13;15 5;21 48;29 39;33 44;15 19;16 19;12 17;50 40;22 53;21 36;
        20 30;26 29;40 20;36 26;62 48;67 41;62 35;65 27;62 24;55 20;35 51;30 50;
        45 42;21 45;36 6;6 25;11 28;26 59;30 60;22 22;27 24;30 20;35 16;54 10;50 15;
        44 13;35 60;40 60;40 66;31 76;47 66;50 70;57 72;55 65;2 38;7 43;9 56;15 56;
        10 70;17 64;55 57;62 57;70 64;64 4;59 5;50 4;60 15;66 14;66 8;43 26];%75 cities d'=549.18 by D B Fogel
    for i=1:75
        for j=1:75
            dislist(i,j)=((city75(i,1)-city75(j,1))^2+(city75(i,2)-city75(j,2))^2)^0.5;
        end
    end
    cityn=city75;
end
end
 
%------------------------------------------------
%��Ӧ�Ⱥ���
function F=CalDist(dislist,s)
 
DistanV=0;
n=size(s,2);
for i=1:(n-1)
    DistanV=DistanV+dislist(s(i),s(i+1));
end
DistanV=DistanV+dislist(s(n),s(1));
F=DistanV;
 
end
 
%------------------------------------------------
%��ͼ
function drawTSP(Clist,x,y,p,index)
CityNum=size(Clist,1);
for i=1:CityNum-1
    plot([Clist(x(i),1),Clist(x(i+1),1)],[Clist(x(i),2),Clist(x(i+1),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
    text(Clist(x(i),1),Clist(x(i),2),['  ',int2str(x(i))]);
    text(Clist(x(i+1),1),Clist(x(i+1),2),['  ',int2str(x(i+1))]);
    hold on;
end
plot([Clist(x(CityNum),1),Clist(x(1),1)],[Clist(x(CityNum),2),Clist(x(1),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
title([num2str(CityNum),'����TSP']);
if index==0&&CityNum~=10
    text(1000,1000,['�� ',int2str(p),' ��','  ��̾���Ϊ ',num2str(y)]);
else
    text(1000,1000,['���������������̾��� ',num2str(y),'�� �ڵ� ',num2str(p),' ���ﵽ']);
end
if CityNum==10
    if index==0
        text(1000,1000,['�� ',int2str(p),' ��','  ��̾���Ϊ ',num2str(y)]);
    else
        text(1000,1000,['���������������̾��� ',num2str(y),'�� �ڵ� ',num2str(p),' ���ﵽ']);
    end
end
hold off;
pause(0.001); 
end
