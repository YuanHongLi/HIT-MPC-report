clc;
clf;
close all;
format long;
A=[-4.59 -0.94;1.52 -4.44];
Bu=[2.29;-0.76];
Bd=[2.30;10.67];
C=[1 0;0 1];
T=0.02;
[Ga,Hu]=c2d(A,Bu,T);
[Gb,Hd]=c2d(A,Bd,T);
p=50;
m=50;
% [x1(k+1);x2(k+1)]=Ga[x1(k);x2(k)]+Hu u(k)+Hd d(k)
simtime = 5.0;
x1 = zeros(simtime/T,2);
%%% ����������ɢϵͳ�Ľ�Ծ��Ӧ�����ڲ����ý�Ծģ�ͣ����Ժ�������û���õ��ⲿ�ֽ��
for i=2:1:simtime/T
    x1(i,:)=Ga*[x1(i-1,1);x1(i-1,2)]+Hu;   
end
for i=2:1:simtime/T
    if x1(i,1)>=x1(end,1)*0.95 
        break
    end
end
N_1 = i;
for i=2:1:simtime/T
    if (x1(i,2)-x1(2,2))>=(x1(end,2)-x1(2,2))*0.95 
        break
    end
end
N_2 = i;
N_u = max(N_1,N_2);
s_u = x1(2:N_u,:);

x1 = zeros(simtime/T,2);
for i=2:1:simtime/T
    x1(i,:)=Ga*[x1(i-1,1);x1(i-1,2)]+Hd;   
end
for i=2:1:simtime/T
    if x1(i,1)>=x1(end,1)*0.95 
        break
    end
end
N_1 = i
for i=2:1:simtime/T
    if (x1(i,2)-x1(2,2))>=(x1(end,2)-x1(2,2))*0.95 
        break
    end
end
N_2 = i
N_d = max(N_1,N_2);
s_d = x1(2:N_d,:);
N=length(s_u);
%%%% ������ɢ״̬�ռ�ģ�͵�Լ��Ԥ�����%%%%
%%%���ɶಽԤ�����%%
Sx=[];
Sd=[];
Su=[];
sum1 = zeros(2,2);
sum2 = zeros(2,1);
sum3 = zeros(2,1,m);
Mss = cell([p,p]);
Mss(:) = {zeros([2,2])};

if p<N
    for i=2:1:p
        Mss{i-1,i} = eye(2,2);
    end
end
Mss{end,end} = eye(2);
Mss = cell2mat(Mss);
for i=1:1:p

    sum1 = sum1 + Ga^(i);
    Sx=[Sx;sum1];
    
    sum2 = sum2 + Ga^(i-1)*Hd;
    Sd=[Sd;sum2];
    
    sum3(:,:,1) = sum3(:,:,1)+Ga^(i-1)*Hu
    
    for j = 2:1:i
            sum3(:,:,j) = sum3(:,:,j-1)-(Ga^(i-j+1)*Hu);
    end
    
    temp=[];
    for j=1:1:m
        temp=[temp sum3(:,:,j)];
    end
    Su=[Su;temp];
end
Fy=eye(2*p,2*p)*0.2;
Fu=eye(m,m)* 3.0;
ymin = [-1;-0.85];
ymax = [ 1; 0.85];
Ymin = cell([p,1]);
Ymax = cell([p,1]);
Ymin(:) = {ymin};
Ymax(:) = {ymax};
Ymin = cell2mat(Ymin);
Ymax = cell2mat(Ymax);
H = (Su')*(Fy')*Fy*Su+(Fu')*Fu;
L = tril(ones(m));
C=[-Su;Su];
yk_=zeros(2,1);
yk = zeros(2*p,1);
x_r_1 = zeros(simtime/T,1);
x_r_2 = zeros(simtime/T,1);
Rp = ones(2*p,1)*0.0;
u = zeros(simtime/T,1);
uk_co = zeros(1,m);
uk_co(1)=1;
deltax = zeros(2,1);
deltay = zeros(2,1);
for i= 1:1:simtime/T
    if i == 1
        dk = 0.1;
    else
        dk = 0.0;
    end
    %Ep = Rp-Mss*yk-Sd*dk;
    if i==1
        Ep = Rp-Sd*dk;
    else
        Ep = Rp-Sx*deltax-yk-Sd*dk;
    end
    
    Gk = 2*Su'*(Fy')*Fy*Ep;
    if i==1
        bk = [Sd*dk-Ymax;Ymin-Sd*dk];
    else
        bk = [Sx*deltax+yk+Sd*dk-Ymax;Ymin-Sx*deltax-yk-Sd*dk];
    end
    
    du = quadprog(H,(-Gk)',-C,-bk);
    if i == 1
    	u(i) = du(1);
    else
        u(i) = u(i-1) + du(1);
    end
    %%% ���ݼ�����Ŀ�����������Ԥ��
    if i==1
         yk = Su*uk_co'*du(1)+Sd*dk;
    else
         yk = Sx*deltax+yk+Su*uk_co'*du(1)+Sd*dk;
    end
    %%% ���ݼ�����Ŀ�����������״̬��������Ϊ������ȫ׼ȷ����L=0
    deltax = Ga*deltax+Hu*du(1)+Hd*dk;
    deltay = eye(2,2)*deltax + deltay;
    x_r_1(i) = deltay(1);
    x_r_2(i) = deltay(2);
    yk(1:2:end) = deltay(1);
    yk(2:2:end) = deltay(2);

end
figure(1)
plot((1:1:simtime/T)*T,u,'-','LineWidth',3);
xlabel('t/s');
title('Control Inpout of u');
legend('u')
hh=findobj('tag','legend') %|
set(hh,'fontsize',10) %| ����legend�ֺŴ�С
set(hh,'LineWidth',1.0); %| ����ͼ���߿�
set(gca,'linewidth',1.5) %| ����ͼ����߿���߿�1.5
set(gca,'box','off') %| ȥͼ�����
figure(2)
plot((1:1:simtime/T)*T,x_r_1, '-','LineWidth',3);
title('\beta');
legend('\beta');
xlabel('t/s');
hh=findobj('tag','legend') %|
set(hh,'fontsize',10) %| ����legend�ֺŴ�С
set(hh,'LineWidth',1.0); %| ����ͼ���߿�
set(gca,'linewidth',1.5) %| ����ͼ����߿���߿�1.5
set(gca,'box','off') %| ȥͼ�����
figure(3)
plot((1:1:simtime/T)*T,x_r_2, '-','LineWidth',3);
title('r');
legend('r');
xlabel('t/s');
hh=findobj('tag','legend') %|
set(hh,'fontsize',10) %| ����legend�ֺŴ�С
set(hh,'LineWidth',1.0); %| ����ͼ���߿�
set(gca,'linewidth',1.5) %| ����ͼ����߿���߿�1.5
set(gca,'box','off') %| ȥͼ�����
