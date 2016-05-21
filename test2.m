%%时序分析
clear all;clc;close all;

d=2;T=0.02;
s=[0 0 1*ones(1,3) 2*ones(1,3) 0*ones(1,6) 3 0*ones(1,5) 3 1*ones(1,6) 2*ones(1,4) 0*ones(1,10) 3 1*ones(1,8) 2*ones(1,5) 0*ones(1,10) 3 1*ones(1,5) 2*ones(1,5) 1*ones(1,8)...
   0*ones(1,8) 2*ones(1,5) 0*ones(1,8) 1*ones(1,5) 3 2*ones(1,5) 0*ones(1,5)];
% s=[0 0 0*ones(1,3) 2*ones(1,3) 3 1*ones(1,4) 3 3 3 3 0*ones(1,5) 3 3 3 1 1 1 0 0]
ll=length(s)
st(1)=0;
s1=s;s2=s;s3=s;

%%时延和时序乱序情况下的控制信号选取一
% for i=0:ll-2
%     if s(i+2)==0
%         st(i+2)=0;
%     elseif s(i+1)<2 & s(i+2)>0
%         st(i+2)=1;
%     else
%         st(i+2)=2;
%     end
%     if st(i+1)==st(i+2)
%         l=l;
%     else
%         l=l+1;
%     end
% end

%%时延和时序乱序情况下的控制信号选取二
% tic
% for i=3:ll
%     for j=0:d
%         if s(i-j)-j<=0
%             sd(i)=j;
%             break
%         end
%     end
% end
% toc

%%时延和时序乱序以及丢包情况下的控制信号选取
for i=1:ll
    if s(i)<3
        p(i)=0;
        pd=0;
    else
        pd=pd+1; %计数器
        p(i)=pd; %连续丢包数序列
        s(i)=2; %时延序列
    end
end

for i=3:ll
    for j=0:d
        if s(i-j)-(j+p(i-j))<=0
            sd(i)=j+p(i-j);
            break
        end
    end
end

m=1;
for i=1:ll
    if s1(i)>2
        ss(m)=i;
        s1(i)=2;
        m=m+1;
    end
end
mm=1;
for i=1:ll-1
    if s3(i)-s3(i+1)>1 & s3(i)<=d
        ss3(mm)=i;
        mm=mm+1;
    end
end

t=0:ll-1;
figure
plot(t,s2*T,'k.')
axis([0 ll-1 -0.5*T 2.5*T])
hold on
plot(ss-1,s1(ss)*T,'r>')
hold on
plot(ss3,s3(ss3+1)*T,'b*')
legend('time delay','packet dropout','packet disordering')
xlabel('Time step k','Fontsize',12)
ylabel('$\tau_k(s)$','Fontsize',12,'interpreter','latex')
grid on

l=0;ls=0;
for i=0:ll-2
    if sd(i+1)~=sd(i+2)
        l=l+1;
    end
    if s(i+1)~=s(i+2)
        ls=ls+1;
    end
end
l
ls

Ap=[0 1;31.5397 0];Bp=[0 -2.6708]';
[A,B]=c2d(Ap,Bp,T);
X(:,1)=[0.1 0 0 0 0];
Xe(1)=norm(X(:,1));
% K = [23.4804    4.1934   -0.1591   -0.1383   -0.1292];%T=0.03
K = [22.5181    4.0244   -0.0948   -0.0835   -0.0888];

[m1,n1]=size(A);
[m2,n2]=size(B);

G0=[A            zeros(m1,n2) zeros(m1,n2) zeros(m1,n2)
    zeros(n2,n1) zeros(n2,n2) zeros(n2,n2) zeros(n2,n2)
    zeros(n2,n1) eye(n2)      zeros(n2,n2) zeros(n2,n2)
    zeros(n2,n1) zeros(n2,n2) eye(n2)      zeros(n2,n2)];
G1=[A            B            zeros(m1,n2) zeros(m1,n2)
    zeros(n2,n1) zeros(n2,n2) zeros(n2,n2) zeros(n2,n2)
    zeros(n2,n1) eye(n2)      zeros(n2,n2) zeros(n2,n2)
    zeros(n2,n1) zeros(n2,n2) eye(n2)      zeros(n2,n2)];
G2=[A            zeros(m1,n2) B            zeros(m1,n2)
    zeros(n2,n1) zeros(n2,n2) zeros(n2,n2) zeros(n2,n2)
    zeros(n2,n1) eye(n2)      zeros(n2,n2) zeros(n2,n2)
    zeros(n2,n1) zeros(n2,n2) eye(n2)      zeros(n2,n2)];
G3=[A            zeros(m1,n2) zeros(m1,n2)            B
    zeros(n2,n1) zeros(n2,n2) zeros(n2,n2) zeros(n2,n2)
    zeros(n2,n1) eye(n2)      zeros(n2,n2) zeros(n2,n2)
    zeros(n2,n1) zeros(n2,n2) eye(n2)      zeros(n2,n2)];
H0=[B;eye(n2);zeros(n2,n2);zeros(n2,n2)];
H1=[zeros(m2,n2);eye(n2);zeros(n2,n2);zeros(n2,n2)];
H2=[zeros(m2,n2);eye(n2);zeros(n2,n2);zeros(n2,n2)];
H3=[zeros(m2,n2);eye(n2);zeros(n2,n2);zeros(n2,n2)];

emax=1294;emin=0.3511;
for i=1:ll
    if sd(i)==0
        X(:,i+1)=(G0+H0*K)*X(:,i);
        u(i)=K*X(:,i);
    elseif sd(i)==1
        X(:,i+1)=(G1+H1*K)*X(:,i);
        u(i)=K*X(:,i);
    elseif sd(i)==2
        X(:,i+1)=(G2+H2*K)*X(:,i);
        u(i)=K*X(:,i);
    else
        X(:,i+1)=(G3+H3*K)*X(:,i);
        u(i)=K*X(:,i);
    end
    Xe(i+1)=sqrt(1294/0.3511)*0.992^i*Xe(1);
    XX(i+1)=norm(X(:,i+1));
end

t=0:ll-1;
figure
subplot(2,1,1)
plot(t,X(1,1:ll),'k-',t,X(2,1:ll),'r-.','Linewidth',1.5)
legend('x1','x2')
xlabel('Time step k','Fontsize',12)
ylabel('State trajectories','Fontsize',12)
axis([0 ll-1 -0.18 0.12])
% figure
subplot(2,1,2)
% plot(t,u(1:ll),'k-')
stairs(t,u,'k-','Linewidth',1.5)
xlabel('Time step k','Fontsize',12)
ylabel('Control inputs','Fontsize',12)
axis([0 ll-1 0 2.5])

figure
plot(t,Xe(2:ll+1),'r',t,XX(2:ll+1),'k')


