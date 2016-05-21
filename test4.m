%%乱序比较
clear all;clc;close all;

d=2;T=0.02;
% s=[0 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 0 0*ones(1,50)];
s=[0 0 0 rand(1,998)];
ll=length(s);
for i=1:ll
    if s(i)>0.5
        s(i)=2;
    elseif s(i)<0.4
        s(i)=0;
    else
        s(i)=1;
    end
end
st(1)=0;
s1=s;s2=s;s3=s;

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

X1(:,1)=[0.1 0]';Kk=[30 1];
X(:,1)=X1(:,1);
for i=1:ll
    X1(:,i+1)=A*X1(:,i)+B*Kk*X1(:,i-s(i));
    X(:,i+1)=A*X(:,i)+B*Kk*X(:,i-sd(i));
end

t=0:ll-1;
figure
subplot(2,1,1)
plot(t,X(1,1:ll),'k-',t,X1(1,1:ll),'r-.','Linewidth',1.5)
legend('without packet disordering','with packet disordering')
xlabel('Time step k','Fontsize',12)
ylabel('x1','Fontsize',12)
axis([0 ll-1 -0.15 0.3])
subplot(2,1,2)
plot(t,X(2,1:ll),'k-',t,X1(2,1:ll),'r-.','Linewidth',1.5)
legend('without packet disordering','with packet disordering')
xlabel('Time step k','Fontsize',12)
ylabel('x2','Fontsize',12)
axis([0 ll-1 -1 1.7])
% 639 550
