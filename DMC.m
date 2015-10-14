clear all;clc;close all;
%����DMC����
A_model=[-4 -0.03;0.75 -10];
B_model=[2 0]';
C_model=[0 1];

Ts=0.04;%����ʱ��
N=40;%ģ�ͳ���
P=15;%Ԥ��ʱ��
M=5;%����ʱ��
control_R=0.001;%����Ȩ����ϵ��
error_Q=1;%���Ȩ����ϵ��
correction_h=1;%У��ϵ��
Sv=1;%�趨ֵ���ο��켣

t=[0:Ts:N*Ts];
y0=step(ss(A_model,B_model,C_model,0),t);

%���ƾ���A��A0,���ƾ���K���¶γ���ó���
A=zeros(P,M);
a=zeros(N,1);
for i=1:N
    a(i)=y0(i+1);%y0Ϊϵͳ��Ծ��Ӧ��ϵ��
end
for i=1:P
    for j=1:M
        if i-j+1>0
            A(i,j)=a(i-j+1);
        end
    end
end

K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
y=zeros(N,1);%���ֵ
u=zeros(N,1);%�����������ʼ��
e=zeros(N,1);%�������ʼ��
A0=zeros(P,N-1);
%A0Ϊϵͳģ�ͣ����ڼ���׼ȷ�����ֵ
for i=1:P
    for j=(N-2):-1:1
        if (N-j+1+i-1)<=N
            A0(i,j)=a(N-j+1+i-1)-a(N-j+i-1);
        else
            A0(i,j)=0;
        end
    end
    A0(i,N-1)=a(i+1);
end

%ѭ�����֣�������ʼֵΪ0
%DMC����
for k=2:N
    Uk_1=zeros(N-1,1);
    for i=1:(N-1)
        if k-N+i<=0
            Uk_1(i)=0;
        else
            Uk_1(i)=u(k-N+i);
        end
    end
    Y0=A0*Uk_1;
    e(k)=y(k-1)-Y0(1);%�������
    Ysk=zeros(P,1);
    for i=1:P
        Ysk(i)=Sv;%Ysk=[Sv Sv ... Sv]'
    end
    Ek=zeros(P,1);
    for i=1:P
        Ek(i)=e(k);%Ek=[e(k) e(k) ...e(k)]'
    end
    dersu=K*(Ysk-Y0-Ek);%�����������
    for i=1:M
        if k+i-1<=N
            u(k+i-1)=u(k+i-1-1)+dersu(i);%���������£������Ż�
        end
    end
    temp=0;
    for j=1:N-1
        if k-j<=0
            temp;
        else
            if k-j-1<=0  
                temp=temp+a(j)*u(k-j);
            else
                temp=temp+a(j)*(u(k-j)-u(k-j-1));%Ԥ�����
            end
        end
    end
    if k-N<=0
        y(k)=temp+e(N)*correction_h;%��������Ԥ�����У��
    else
        y(k)=temp+a(N)*u(k-N)+e(N)*correction_h;%��������Ԥ�����У��
    end
end

%��ͼ��ʾ���
t=Ts.*(1:N);
subplot(2,1,1);
plot(t,y,'.-');
legend('y','Location','Best');
title('�������');
xlabel('Time')
ylabel('���')
grid on
subplot(2,1,2);
plot(t,u,'.-');
legend('��������u','Location','Best');
title('��������');
xlabel('Time')
ylabel('���')
grid on
hold on
    
    