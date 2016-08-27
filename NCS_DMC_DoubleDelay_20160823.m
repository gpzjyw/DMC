%  2016-08-23：现在做如下改进：1.考虑双边时延的情况，且都可能出现乱序的问题；
%                             2.乱序的处理，执行器选取最新控制信号，控制器也选用最新的输出值；
%                             3.控制器的校正部分，所用的控制增量根据时延 传感器-控制器 时延确定；
%                             4.控制器采用事件驱动，收到一个输出值驱动一次事件，但选用的输出值为最新的输出值；
% ********************************************************************************** %

clear all;clc;close all;

d=4; % 最大时延为 d*T
T=0.04; % 采样周期
totalStep=50;% 仿真总步长
timeSequence = d*T*[0 rand(1,totalStep-d-1)-0.2 zeros(1,d)]; % 每个时刻数据包的时延大小（发生时延的概率为80%）
timeSequenceInt = timeSequence; % 缓冲区的存在，将时延转化为固定周期数的时延：0,T,2T,...,dT
timeSequenceLength = length(timeSequence);
controlSequence=zeros(2*(d+1),totalStep);

% 将时延序列转义为0,1,2,3....
for i=1:totalStep
    if timeSequence(i)<=0
        timeSequence(i)=0;
        timeSequenceInt(i)=0;
    end
    for j=1:d
        if timeSequence(i)<=j*T
            if timeSequence(i)>(j-1)*T
                timeSequenceInt(i)=j;
                break
            end
        end
    end
end

% 设定丢包的概率，将发生丢包的时刻存入dropoutPoint，并将丢包时刻的时延设置为(d+1)T，发生丢包的概率为10%
dropoutRate=[0 rand(1,totalStep-d-1) zeros(1,d)];
dropoutPoint=[]; % 丢包时刻的集合
for i=1:totalStep
    if dropoutRate(i)>=0.9
        timeSequenceInt(i)=d+1;
        dropoutPoint=[dropoutPoint,i];
    end
end

% 确定每个时刻之前的连续丢包数
continuousDropoutNum=zeros(1,totalStep);
for i=1:length(dropoutPoint)
    if continuousDropoutNum(dropoutPoint(i))~=0
        continuousDropoutNum(dropoutPoint(i)+1)=continuousDropoutNum(dropoutPoint(i))+1;
    else
        continuousDropoutNum(dropoutPoint(i)+1)=1;
    end
end

% 根据时延确定乱序时刻的集合
disorderingPoint=[]; % 乱序时刻的集合
for k=1:totalStep
    
    if k<=d+1
        dValue=k-1;
    end
    
    if timeSequenceInt(k)~=0
        if k>3
            for i=1:(dValue-1)
                for j=(i+1):dValue
                    if (timeSequenceInt(k-j)-timeSequenceInt(k-i))>=2
                        if timeSequenceInt(k-j)==j
                            flag=1;
                        end
                    end
                end
            end
        end
    end
    
	if flag==1
        disorderingPoint=[disorderingPoint,k];
        flag=0;
	end    
    
end

%取出未发生丢包的时刻的时延，存入timeSequenceIntDraw中；并将未丢包的时刻存入delayPoint中
timeSequenceIntDraw=[];
j=1;
delayPoint=[]; % 时延时刻的集合
for i=1:totalStep
    if i==dropoutPoint(j)
        j=j+1;
        if j>=length(dropoutPoint)
            j=length(dropoutPoint);
        end
    else
        timeSequenceIntDraw=[timeSequenceIntDraw,timeSequenceInt(i)];
        delayPoint=[delayPoint,i];
    end
end

% 计算sigma(k)，执行器根据该值选取最新的控制量
sigma=zeros(1,totalStep);
for k=1:totalStep
    if k<=d+1
        dValue=k;
    end
    for i=0:dValue
        if timeSequenceInt(k-i)<=i
            sigma(k)=i;
            break;
        end
    end
end

% 绘图
figure(1)
t=T*(1:totalStep);
subplot(2,1,1);
plot(delayPoint*T,timeSequenceIntDraw*T,'*k',disorderingPoint*T,zeros(length(disorderingPoint),1),'ok',dropoutPoint*T,d*T,'xk');
axis([0 totalStep*T -0.5*T (d+1)*T]);
legend('时延','发生乱序的时刻','丢包的数据');
xlabel('时间k');
ylabel('\tau _k');
grid on
hold on
subplot(2,1,2);
plot((1:totalStep)*T,sigma,'*k')
xlabel('时间k');
ylabel('\sigma(k)');
grid on
hold on



% ********************************************************************************** %
% 被控对象的状态空间方程
A_model=[-4 -0.03;0.75 -10];
B_model=[2 0]';
C_model=[0 1];

% 预定义DMC计算中所需要的相关参数
N=40; % 模型时域
P=15; % 预测时域
M=5; % 控制时域
control_R=0.001; % 控制权矩阵系数
error_Q=1; % 误差权矩阵系数
correction_h=1; % 校正系数
Sv=1; % 设定值，参考轨迹
simulationTime=0+totalStep;

% 系统的离散状态空间方程
[A_model_discrete,B_model_discrete]=c2d(A_model,B_model,T);
C_model_discrete=C_model;

% 无时延系统的阶跃响应，获得阶跃响应系数
temp_response=dstep(A_model_discrete,B_model_discrete,C_model_discrete,0);
stepResponse(:,1)=[zeros(N,1),eye(N),zeros(N,length(temp_response)-(N+1))]*temp_response;




























