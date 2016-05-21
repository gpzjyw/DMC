clear all;clc;close all;

d=4;%���ʱ��Ϊ d*T
T=0.04;%��������
% timeSequence = d*T*[0 0.6 0 0.2 0 0.7 0 0.1 0 0.8 0 0.3 0 0.9 0 0.1 0 0.5 0 0.1 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0.2 0];%
% timeSequence= d*T*[zeros(1,50)];%��ʱ�ӵ����绯ģ��
timeSequence = d*T*[0 rand(1,47)-0.2 0 0];
timeSequenceInt = timeSequence;
timeSequenceLength = length(timeSequence);
controlSequence=zeros(2*(d+1),timeSequenceLength);

%��ʱ������ת��Ϊ0,1,2,3......
for i=1:timeSequenceLength
    if timeSequence(i)<0
        timeSequence(i)=0;
    end
    for j=0:d
        if timeSequence(i)<=j*T
            if timeSequence(i)>(j-1)*T
                timeSequenceInt(i)=j;
                break
            end
        end
    end
end

%���г�һ��2(d+1)*timeSequenceLength�ľ����� ��1,3,5,...,(2d+1)�е�k�� ��ʾkʱ���Ƿ��յ� k,(k-1),(k-2)...,(k-d)ʱ�̵Ŀ�����
for i=1:timeSequenceLength
    for j=0:d
        if timeSequenceInt(i)==j
            controlSequence(2*j+1,i+j)=1;
        end
    end
end

figure(1)
t=T*(1:timeSequenceLength);
subplot(2,1,1);
% plot(t,timeSequenceInt*T,'o',t,timeSequence,'.');
plot(t,timeSequenceInt*T,'.');
title('�������ʱ�ӷֲ�');
xlabel('Time/s');
ylabel('ʱ��');
grid on
hold on

subplot(2,1,2);
plot(t,timeSequenceInt,'*')
xlabel('Time/s');
ylabel('timeSequenceInt/ʱ�ӵ�������');
grid on
hold on


A_model=[-4 -0.03;0.75 -10];
B_model=[2 0]';
C_model=[0 1];
N=40;%ģ��ʱ��
P=15;%Ԥ��ʱ��
M=5;%����ʱ��
control_R=0.001;%����Ȩ����ϵ��
error_Q=1;%���Ȩ����ϵ��
correction_h=1;%У��ϵ��
Sv=1;%�趨ֵ���ο��켣
simulationTime=0+timeSequenceLength;

%��ý�Ծ��Ӧ��ϵ��
t=[0:T:N*T];
stepResponse=step(ss(A_model,B_model,C_model,0),t);

% figure(3)
% plot(t,stepResponse);
% hold on

A=zeros(P,M);
a=zeros(N,1);
for i=1:N
	a(i)=stepResponse(i+1);%aΪϵͳ��Ծ��Ӧ��ϵ��
end
for i=1:P
	for j=1:M
        if i-j+1>0 
            A(i,j)=a(i-j+1);
        end
    end
end
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


%����ѡȡk+1ʱ�̵����¿�����������ʱ���޿��������������һʱ�̿���������reordering
for k=1:timeSequenceLength-1 %NӦ����timeSequenceLength����һ��
    
    if k==1
        Y_setValue=ones(P,1)*Sv;%Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        controlValue=zeros(1,simulationTime);%��������ʼ��
        actualControlValue=controlValue;%ʵ�ʿ�������ʼ��
        controlValue(1)=1;
        actualControlValue=controlValue(1);
        
        controlIncrement=zeros(M,simulationTime);%����������ʼ��
        controlIncrement(:,1)=zeros(M,1);
        
        h=ones(N,1);%У������
        
        error=zeros(1,simulationTime);
        
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=0;
        
             
        Y_predictedValue=zeros(N,simulationTime);%Y��Ԥ��ֵ��ʼ��
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);%����ĳ�ʼԤ��ֵ
        
    end
    
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k));%�����������
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1);%���ݿ����������������
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+a*controlIncrement(1,k+1);%���ݿ�����������(k+1��ʱ��y��Ԥ��ֵ 

        if k<=d
            dValue=k;
        end
        %����ʱ�ӣ����������洢��controlSequence��Ӧ��λ��
        for temp1=0:dValue
           if timeSequenceInt(k+1)==temp1
               controlSequence(2*(temp1+1),k+1+temp1)=controlValue(k+1);
           end
        end
        
%         for temp1=0:dValue
%            if timeSequenceInt(k+1)==temp1
%                controlSequence(2*(temp1+1),k+1+temp1)=controlIncrement(1,k+1);
%            end
%         end
        
        actualSingal=0;
        for temp2=0:dValue
            if controlSequence(2*temp2+1,k+1)==1
                actualControlValue(k+1)=controlSequence(2*temp2+2,k+1);
%                 actualControlValue(k+1)=actualControlValue(k)+controlSequence(2*temp2+2,k+1);
                actualSingal=1;
                break
            end
        end
        if actualSingal==0
            actualControlValue(k+1)=actualControlValue(k);
        end
        
        Uk_1=zeros(N-1,1);
        for i=1:(N-1)
            if k+1-N+i<0
                Uk_1(i)=0;
            else
                Uk_1(i)=actualControlValue(k+1-N+i+1);
            end
        end
        Y0=A0*Uk_1;
        Y_outputValue(k+1)=Y0(1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1);%����(k+1)ʱ��y��ʵ��ֵ��Ԥ��ֵ������������У��
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1);%��������ʽ���Ԥ�ⷽ��������(k+1)ʱ�̵�Ԥ��ֵ
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1];%������λ����S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

figure(2)
subplot(2,1,1)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-');
title('�������');
xlabel('Time/s');
ylabel('���');
grid on
hold on

subplot(2,1,2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'-');
title('������');
xlabel('Time/s');
ylabel('���');
grid on
hold on


%�����˲���disordering�Ŀ���������reordering
for k=1:timeSequenceLength-1 %NӦ����timeSequenceLength����һ��
    
    if k==1
        Y_setValue=ones(P,1)*Sv;%Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        controlValue=zeros(1,simulationTime);%��������ʼ��
        actualControlValue=controlValue;%ʵ�ʿ�������ʼ��
        controlValue(1)=1;
        actualControlValue=controlValue(1);
        
        controlIncrement=zeros(M,simulationTime);%����������ʼ��
        controlIncrement(:,1)=zeros(M,1);
        
        h=ones(N,1);%У������
        
        error=zeros(1,simulationTime);
        
%         X_stateValue=zeros(2,simulationTime);
%         X_stateValue(:,1)=[0 0]';
        
        Y_outputValue=zeros(1,simulationTime);
%         Y_outputValue(1)=C_discreteModel*X_stateValue(:,1);
        Y_outputValue(1)=0;
        
             
        Y_predictedValue=zeros(N,simulationTime);%Y��Ԥ��ֵ��ʼ��
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);%����ĳ�ʼԤ��ֵ
        
    end
    
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k));%�����������
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1);%���ݿ����������������
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+a*controlIncrement(1,k+1);%���ݿ�����������(k+1��ʱ��y��Ԥ��ֵ

        if k<=d
            dValue=k;
        end
        
        temp3=zeros(dValue+1,1);
        for temp2=0:dValue
            temp3(temp2+1)=timeSequenceInt(k+1-temp2)-temp2;
        end
%         minIndex,��Сֵ����
%         minValue����Сֵ
        for temp2=0:dValue
            if temp3(temp2+1)<=0
                minIndex=temp2;
                minValue=temp3(minIndex+1);
                break
            end
        end
        
        actualControlValue(k+1)=controlValue(k+1-minIndex);
%         actualControlValue(k+1)=actualControlValue(k)+controlIncrement(1,k+1-minIndex);
        
        Uk_1=zeros(N-1,1);
        for i=1:(N-1)
            if k+1-N+i<0
                Uk_1(i)=0;
            else
                Uk_1(i)=actualControlValue(k+1-N+i+1);
            end
        end
        Y0=A0*Uk_1;
        Y_outputValue(k+1)=Y0(1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1);%����(k+1)ʱ��y��ʵ��ֵ��Ԥ��ֵ������������У��
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1);%��������ʽ���Ԥ�ⷽ��������(k+1)ʱ�̵�Ԥ��ֵ
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1];%������λ����S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

figure(2)
subplot(2,1,1)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'.--');
% legend('without reordering','with reordering');
title('�������');
xlabel('Time/s');
ylabel('���');
grid on
hold on

subplot(2,1,2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'.--');
% legend('without reordering','with reordering');
title('������');
xlabel('Time/s');
ylabel('���');
grid on
hold on

for k=1:timeSequenceLength-1 %NӦ����timeSequenceLength����һ��
    
    if k==1
        Y_setValue=ones(P,1)*Sv;%Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        controlValue=zeros(1,simulationTime);%��������ʼ��
        actualControlValue=controlValue;%ʵ�ʿ�������ʼ��
        controlValue(1)=1;
        actualControlValue=controlValue(1);
        
        controlIncrement=zeros(M,simulationTime);%����������ʼ��
        controlIncrement(:,1)=zeros(M,1);
        
        h=ones(N,1);%У������
        
        error=zeros(1,simulationTime);
        
%         X_stateValue=zeros(2,simulationTime);
%         X_stateValue(:,1)=[0 0]';
        
        Y_outputValue=zeros(1,simulationTime);
%         Y_outputValue(1)=C_discreteModel*X_stateValue(:,1);
        Y_outputValue(1)=0;
        
             
        Y_predictedValue=zeros(N,simulationTime);%Y��Ԥ��ֵ��ʼ��
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);%����ĳ�ʼԤ��ֵ
        
    end
    
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k));%�����������
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1);%���ݿ����������������
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+a*controlIncrement(1,k+1);%���ݿ�����������(k+1��ʱ��y��Ԥ��ֵ

        if k<=d
            dValue=k;
        end
        
        temp3=zeros(dValue+1,1);
        for temp2=0:dValue
            temp3(temp2+1)=timeSequenceInt(k+1-temp2)-temp2;
        end
%         minIndex,��Сֵ����
%         minValue����Сֵ
        for temp2=0:dValue
            if temp3(temp2+1)<=0
                minIndex=temp2;
                minValue=temp3(minIndex+1);
                break
            end
        end
        
        actualControlValue(k+1)=controlValue(k)+controlIncrement(1,k+1-minIndex);
        
        Uk_1=zeros(N-1,1);
        for i=1:(N-1)
            if k+1-N+i<0
                Uk_1(i)=0;
            else
                Uk_1(i)=actualControlValue(k+1-N+i+1);
            end
        end
        Y0=A0*Uk_1;
        Y_outputValue(k+1)=Y0(1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1);%����(k+1)ʱ��y��ʵ��ֵ��Ԥ��ֵ������������У��
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1);%��������ʽ���Ԥ�ⷽ��������(k+1)ʱ�̵�Ԥ��ֵ
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1];%������λ����S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

figure(2)
subplot(2,1,1)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'.--');
legend('without reordering','with reordering1','with reordering2');
title('�������');
xlabel('Time/s');
ylabel('���');
grid on
hold on

subplot(2,1,2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'.--');
legend('without reordering','with reordering1','with reordering2');
title('������');
xlabel('Time/s');
ylabel('���');
grid on
hold on