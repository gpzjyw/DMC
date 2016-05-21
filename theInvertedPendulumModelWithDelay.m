clear all;clc;close all;

d=4;%最大时延为 d*T
T=0.04;%采样周期
% timeSequence = d*T*[0 0.6 0 0.2 0 0.7 0 0.1 0 0.8 0 0.3 0 0.9 0 0.1 0 0.5 0 0.1 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0 0.6 0 0.2 0.2 0];%
timeSequence= d*T*[zeros(1,500)];%无时延的网络化模型
% timeSequence = d*T*[0 rand(1,47)-0.2 0 0];
timeSequenceInt = timeSequence;
timeSequenceLength = length(timeSequence);
controlSequence=zeros(2*(d+1),timeSequenceLength);

%将时间序列转义为0,1,2,3......
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

%排列成一个2(d+1)*timeSequenceLength的矩阵，用 第1,3,5,...,(2d+1)行第k列 表示k时刻是否收到 k,(k-1),(k-2)...,(k-d)时刻的控制量
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
title('各个点的时延分布');
xlabel('Time/s');
ylabel('时延');
grid on
hold on

subplot(2,1,2);
plot(t,timeSequenceInt,'*')
xlabel('Time/s');
ylabel('timeSequenceInt/时延的周期数');
grid on
hold on


A_model=[0 1 0 0;0 0 0 0;0 0 0 1;0 0 29.4 0];
B_model=[0 1 0 3]';
C_model=[1 0 0 0;0 0 1 0];

N=40;%模型时域
P=15;%预测时域
M=5;%控制时域
control_R=0.001;%控制权矩阵系数
error_Q=1;%误差权矩阵系数
correction_h=1;%校正系数
Sv=0;%设定值，参考轨迹
simulationTime=0+timeSequenceLength;

%获得阶跃响应的系数
t=[0:T:N*T];
stepResponse=step(ss(A_model,B_model,C_model,0),t);

figure(3)
plot(t,stepResponse);
hold on

%{
A=zeros(P,M);
a=zeros(N,1);
for i=1:N
	a(i)=stepResponse(i+1);%a为系统阶跃响应的系数
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

%处理了产生disordering的控制量，有reordering
for k=1:timeSequenceLength-1 %N应该与timeSequenceLength保持一致
    
    if k==1
        Y_setValue=ones(P,1)*Sv;%Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        controlValue=zeros(1,simulationTime);%控制量初始化
        actualControlValue=controlValue;%实际控制量初始化
        controlValue(1)=1;
        actualControlValue=controlValue(1);
        
        controlIncrement=zeros(M,simulationTime);%控制增量初始化
        controlIncrement(:,1)=zeros(M,1);
        
        h=ones(N,1);%校正向量
        
        error=zeros(1,simulationTime);
        
%         X_stateValue=zeros(2,simulationTime);
%         X_stateValue(:,1)=[0 0]';
        
        Y_outputValue=zeros(1,simulationTime);
%         Y_outputValue(1)=C_discreteModel*X_stateValue(:,1);
        Y_outputValue(1)=0;
        
             
        Y_predictedValue=zeros(N,simulationTime);%Y的预测值初始化
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);%输出的初始预测值
        
    end
    
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k));%计算控制增量
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1);%根据控制增量计算控制量
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+a*controlIncrement(1,k+1);%根据控制增量计算(k+1）时刻y的预测值

        if k<=d
            dValue=k;
        end
        
        temp3=zeros(dValue+1,1);
        for temp2=0:dValue
            temp3(temp2+1)=timeSequenceInt(k+1-temp2)-temp2;
        end
%         minIndex,最小值索引
%         minValue，最小值
        for temp2=0:dValue
            if temp3(temp2+1)<=0
                minIndex=temp2;
                minValue=temp3(minIndex+1);
                break
            end
        end
        
        actualControlValue(k+1)=controlValue(k+1-minIndex);
        
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
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1);%根据(k+1)时刻y的实际值和预测值计算误差，以用于校正
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1);%运用启发式误差预测方法，修正(k+1)时刻的预测值
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1];%给出移位矩阵S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

figure(2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'-');
legend('without reordering','with reordering');
title('控制量');
xlabel('Time/s');
ylabel('振幅');
grid on
hold on

figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-');
legend('without reordering','with reordering');
title('输出曲线');
xlabel('Time/s');
ylabel('振幅');
grid on
hold on

%}