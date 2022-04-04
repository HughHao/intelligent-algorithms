function [G_k,G_best_route,G_best_length,K,best_route,best_length,length_ave]=VRP_original(C,NC_max,m,Alpha,Beta,Rho,Q,W,t) %#ok<INUSD,STOUT>
%%
%--------------------------------------------------------------------------
%%����˵��
%%G_K:�������ų�����Ŀ
%%G_best_route����������·��
%%G_best_length����������·�߳���
%%K�����ų�����Ŀ
%%best_route������·��
%%best_length������·�߳���
%%length_ave������ƽ������
%%C��DC�͹���������
%%NC_max������������
%%m��������Ŀ
%%Alpha����Ҫ��ϵ��
%%Beta���ܼ���ϵ��
%%Rho���ӷ���ϵ��
%%Q����Ϣ���²���
%%W:����������
%%
%%��һ�� ��ʼ��������ϵ��

m=10;Alpha=2;Beta=1;gama=1;Rho=0.1;NC_max=100;Q=15;W=9;qq=0.05; % #ok<NASGU>
C=[0 0
    0 -1
    0 3
    -2 -2
    -3 -3
    3 -1
    -4 0
    -4 -1
    1 -2
    1 -1
    1 3
    3 4
    -3 0
    2 0
    1 -3
    2 -1
    2 1
    1 -4
    -3 2
    -1 -1]; %  #ok<NASGU>
t=[0 1.5 1.8 2 0.8 1.5 1.0 2.5 3.0 1.7 0.6 0.2 2.4 1.9 2.0 0.7 0.5 2.2 3.1 0.1];

%����ֿ�͹�����֮��ľ������

n=size(C,1);%n��ʾ����Ĺ�ģ�����и����� 
D=zeros(n,n);%D��ʾ��ȫͼ�ĸ�Ȩ�ڽӾ��� 
for i=1:n 
for j=1:n 
if i~=j 
D(i,j)=((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5; 
else 
D(i,j)=eps; 
end 
D(j,i)=D(i,j); 
end 
end 

%�����ʡ������
U=zeros(n,n);%U��ʾ����֮������Ӻ�����ֿ�֮�������ܽ�ʡ�ľ���
for i=1:n 
for j=1:n 
if i~=j 
U(i,j)=D(i,1)+D(j,1)-D(i,j);
else 
U(i,j)=eps; 
end 
U(j,i)=U(i,j); 
end 
end 
load_w=0;
Eta=1./D;%EtaΪ�������ӣ�������Ϊ����ĵ��� 
Tau=ones(n,n);%TauΪ��Ϣ�ؾ��� 
Tabu=zeros(m,n);%�洢����¼·�������� 
NC=1;%���������� 
G_best_route=[NC_max,n];%�������·�� 
G_best_length=inf.*ones(NC_max,1);%�������·�ߵĳ��� 
length_ave=zeros(NC_max,1);%����·�ߵ�ƽ������

%%�ڶ����������Ϸŵ�DC��
while NC<=NC_max%ֹͣ����֮һ���ﵽ���������� 
    Tabu(:,1)=randi([1,1],m,1);
%%��������mֻ���ϰ���Ҫ��ķ���ѡ�񹤳������������ 
for i=1:m
    visited=Tabu(i,:);
    visited=visited(visited>0);
    to_visit=setdiff(1:n,visited);
    c_temp=length(to_visit);
    j=1;
       while j<=n
          if ~isempty(to_visit)
    %visit_thisant=[];
    %while to_visit~=[1]
    %for ee=1:(2*c_temp) 
   %while c_temp>=2
   %if length(to_visit)>1
    %visited=Tabu(i,:);
    %to_visit=setdiff(1:n,visited);
    %to_visit=[1,to_visit];
    %c_temp=length(to_visit);
    %visit_thisant=[];
   % if c_temp~=1
        %for b=1:c_temp
            %if (load_w+t(to_visit(b)))<=W
                %visit_thisant=[visit_thisant,to_visit(b)];
            %end
        %end
        %jj=visited(visited>0);
%% ���չ���ѡ��һ�����������ǻص��ֿ�
%for k=1:length(visit_thisant)
    for k=1:length(to_visit)
  x(k)=(Tau(visited(end),to_visit(k))^Alpha)*(Eta(visited(end),to_visit(k))^Beta)*(U(visited(end),to_visit(k))^gama);
    end   
%         ww=rand;
% if ww<qq
%     Select=find(max(x));
%     %Tabu(i,length(visited(visited>0))+1)=to_visit(Select(1)); 
% else
x=x/(sum(x)); 
%������ԭ��ѡȡ��һ������ 
xcum=cumsum(x); 
Select=find(xcum>=rand);
%Tabu(i,length(visited(visited>0))+1)=to_visi(Select(1)); 
% end
if isempty(Select)
    Select=1;
    load_w=load_w+t(Select);
else
load_w=load_w+t(to_visit(Select(1)));
end
if load_w>W
    Select=1;
       j=j-1;
    load_w=0;
    Tabu(i,length(visited)+1)=Select(1);
else
Tabu(i,length(visited)+1)=to_visit(Select(1)); 
end
%if Select(1)==1

%else
    %load_w=load_w+t(Select(1));
%end
          end
    visited=Tabu(i,:);
    visited=visited(visited>0);
    to_visit=setdiff(1:n,visited);
    x=[];
            if visited(end)~=1
   Tabu(i,1:(length(visited)+1))=[visited,1];
            end
        j=j+1;
       end
    load_w=0;
end
%% ���Ĳ���¼�������ֲ���
L=zeros(m,1); 
for i=1:m 
MM=Tabu(i,:); 
R=MM(MM>0);
for j=1:(length(R)-1)
L(i)=L(i)+D(R(j),R(j+1)); 
end 
end 

G_best_length(NC)=min(L); 
pos=find(L==G_best_length(NC)); 
G_best_route(NC,1:length(Tabu(pos(1),:)))=Tabu(pos(1),:);
%%Ӧ��2-opt���������Ž���и���
select=find(G_best_route(NC,:)==1);
FF=[];
LM=0;
for a=1:(length(select)-1)
   y_G_best_route=G_best_route(NC,select(a):select(a+1));
   al=length(y_G_best_route);
   T=zeros((length(select)-1),1);
   for d=1:(al-1)
       T(a)=T(a)+D(y_G_best_route(d),y_G_best_route(d+1));
   end
   for b=2:(al-1)
       for c=(b+1):(al-1)
           DD=y_G_best_route;
           temp1=DD(b);
           temp2=DD(c);
           DD(b)=temp2;
           DD(c)=temp1;
           TT=zeros(1);
           for d=1:(al-1)
       TT=TT+D(DD(d),DD(d+1));
           end
           if TT<T(a)
               T(a)=TT;
            y_G_best_route=DD;
           end
       end
   end
   if a>=2
       y_G_best_route=y_G_best_route(2:al);
   end
   FF=[FF,y_G_best_route];
   LM=LM+T(a);
end
   G_best_length(NC)=LM;
   G_best_route(NC,1:length(FF))=FF;
   FF=[];
   LM=0;
   %%2-opt�Ż����
length_ave(NC)=mean(L); 
NC=NC+1;
%% ���岽������Ϣ��
Delta_Tau=zeros(n,n); 
for i=1:m 
MM=Tabu(i,:); 
R=MM(MM>0);
for j=1:(length(R)-1) 
Delta_Tau(R(j),R(j+1))=Delta_Tau(R(j),R(j+1))+Q/L(i); 
end 
end 
Tau=(1-Rho).*Tau+Delta_Tau;
%% �����������ɱ����� 
Tabu=zeros(m,n); 
load_w=0;
end
%%���߲��������� 
Pos=find(G_best_length==min(G_best_length)); 
best_route=G_best_route(Pos(1),:);
best_length=G_best_length(Pos(1));
best_route=best_route(best_route>0);
best_route
best_length
%% �ڰ˲�������ɢ��ͼ��Ѳ�ι���ͼ

figure(1)
plot([C(best_route,1)],[C(best_route,2)],'-*')
 figure(2)
plot(G_best_length) 
hold on 
plot(length_ave) 

end
