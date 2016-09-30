%  2016-08-23：现在做如下改进：1.考虑双边时延的情况，且都可能出现乱序的问题；
%                             2.乱序的处理，执行器选取最新控制信号，控制器也选用最新的输出值；
%                             3.控制器的校正部分，所用的控制增量根据时延 传感器-控制器 时延确定；
%                             4.控制器采用事件驱动，收到一个输出值驱动一次事件，但选用的输出值为最新的输出值；
%  2016-09-04：现在做如下考虑: 1.去除丢包因素，仅仅考虑时延和乱序的情况；
%                             2.控制器和执行器设置缓冲区，将不确定的时延转化为若干固定时延 0,1,2,3,...,d；
%  2016-09-20：现在做如下考虑: 1.控制器采用事件驱动，执行器采用时间驱动的工作方式；
%                             2.控制器对乱序的数据包进行主动丢弃；
%  2016-09-28：先增加一对照组：1.处理乱序问题，但不对时延进行补偿
%                             2.更改为传输控制量
% ********************************************************************************** %

clear all;clc;close all;

d=4; % 最大时延为 d*T
T=0.04; % 采样周期
totalStep=75;% 仿真总步长
timeSequence = d*T*[0 rand(1,totalStep-d-1)-0.2 zeros(1,d)]; % 每个时刻数据包的时延大小（发生时延的概率为80%）
timeSequenceInt = timeSequence; % 缓冲区的存在，将时延转化为固定周期数的时延：0,T,2T,...,dT

% 总时延
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

% 计算sigma(k)，执行器根据该值选取最新的控制量，适用于：控制器不处理乱序，执行器处理乱序
sigma_SC_noreorder=zeros(1,totalStep);
for k=1:totalStep
    for i=0:(d+1)
        if k-i<=0
            break;
        end
        if timeSequenceInt(k-i)<=i
            sigma_SC_noreorder(k)=i;
            break;
        end
    end
end

% 根据时延，计算任一时刻最新到达的控制量，适用于：控制器不处理乱序，执行器不处理乱序
latestArrive_SC_noreorder=zeros(1,totalStep);
for k=1:totalStep
    flag=0;
    for i=0:d
        if k-i<=0
            break;
        end
        if timeSequenceInt(k-i)==i
           latestArrive_SC_noreorder(k)=i;
           flag=1;
           break;
        end
    end
    if flag==0
        latestArrive_SC_noreorder(k)=latestArrive_SC_noreorder(k-1)+1;
    end
end

% 传感器-控制器时延
% 以 传感器-控制器时延 小于等于 全局时延 为约束，生成 传感器-控制器时延
timeSequence_S2C=zeros(1,totalStep);
timeSequenceInt_S2C=zeros(1,totalStep);
for i=1:totalStep
    timeSequence_S2C(i) = getRandom(0, timeSequence(i));
    timeSequenceInt_S2C(i) = ceil(timeSequence_S2C(i)/T);
end

% 确定 传感器-控制器 乱序的数据包，定义该时刻数据包为丢包（用时延d+1表示）
disorderingPoint_S2C=[]; % 乱序时刻的集合
for k=1:totalStep
    for j=1:d
        if k+j>totalStep
            continue;
        end
        if timeSequence_S2C(k)-timeSequence_S2C(k+j)>j*T
            timeSequenceInt(k)=d+1;
            disorderingPoint_S2C=[disorderingPoint_S2C,k];
            break;
        end
    end
end

% 根据具体到达时间进行重排序
arriveTime_S2C=[[1:totalStep] ; [1:totalStep]*T+timeSequence_S2C];
arriveTimeOrder_S2C=sortrows(arriveTime_S2C', 2)';
% 控制器接收的数据包时戳顺序
timeStampList_S2C=arriveTimeOrder_S2C(1,:);
% 每个采样时刻到达数据包的时间戳，首行向下按顺序排列
timeStampSquare_S2C=zeros(d+1,totalStep);

% 根据时延确定乱序数据包的集合
disorderingPoint=[]; % 乱序时刻的集合
for k=1:totalStep
    
    if k<=d+1
        dValue=k-1;
    end
    
    flag=0;
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
    end
    
end

% 未丢包的数据包的时延，存入timeSequenceIntDraw
% 未丢包的数据包的时间戳，存入delayPoint
% 丢包的数据包的时间戳，存入dropoutPoint
timeSequenceIntDraw=[];
delayPoint=[]; % 时延时刻的集合
dropoutPoint=[]; % 丢包时刻集合
for i=1:totalStep
    if timeSequenceInt(i)==(d+1)
        dropoutPoint=[dropoutPoint,i];
        continue;
    end
    timeSequenceIntDraw=[timeSequenceIntDraw,timeSequenceInt(i)];
    delayPoint=[delayPoint,i];
end

% 计算sigma(k)，执行器根据该值选取最新的控制增量，适用于：控制器处理乱序，执行器处理乱序
sigma_SC_reorder=zeros(1,totalStep);
for k=1:totalStep
    for i=0:(d+1)
        if k-i<=0
            break;
        end
        if timeSequenceInt(k-i)<=i
            sigma_SC_reorder(k)=i;
            break;
        end
    end
end

% 根据时延，计算任一时刻最新到达的控制增量，适用于：控制器处理乱序，执行器不处理乱序
latestArrive_SC_reorder=zeros(1,totalStep);
for k=1:totalStep
    for i=0:d
        if k-i<=0
            break;
        end
        if timeSequenceInt(k-i)==i && timeSequenceInt(k-i)~=(d+1)
           latestArrive_SC_reorder(k)=i;
           break;
        end
    end
end

% 绘图
figure(1)
set(gcf,'Position',[100,400,700,500]);
t=T*(1:totalStep);
subplot(2,1,1);
plot(delayPoint*T,timeSequenceIntDraw*T,'*k',disorderingPoint*T,zeros(length(disorderingPoint),1),'ok',dropoutPoint*T,d*T,'xk');
axis([0 totalStep*T -0.5*T (d+1)*T]);
legend('数据包的时延','执行器发生乱序的时刻','主动丢弃的数据包','Orientation','horizontal','Location','best');
xlabel('时间k');
ylabel('\tau _k');
grid on
hold on
subplot(2,1,2);
plot((1:totalStep)*T,sigma_SC_reorder,'*k')
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

% 时延为1,2,3,...,d时候的阶跃响应系数
% stepResponse 中每列存储着时延为0,1,2,...,2d情况下的阶跃响应系数
for delay=1:(2*d)
    for i=1:N
        accumulation=zeros(2,2);%此处需要根据A_model、B_model的维数确定
        if i<=delay
            for j=0:(i-1)
                accumulation = accumulation + A_model_discrete^j;
            end
        else
            for j=(i-delay):(i-1)
                accumulation = accumulation + A_model_discrete^j;
            end
        end
        stepResponse(i,delay+1) = stepResponse(i,1) - C_model_discrete * accumulation * B_model_discrete;
    end
end

% A_gather中依次排列着时延为0到2d下的动态矩阵
A_gather=zeros(P,M*(2*d+1));
for delay=0:(2*d)
    % 将阶跃响应系数赋给行向量a
    a=zeros(N,1);
    a=stepResponse(:,delay+1); % a为系统阶跃响应的系数
    % 计算不同时延情况下的动态矩阵
    for i=1:P
        for j=1:M
            if i-j+1>0
                A_gather(i,delay*M+j)=a(i-j+1);
            end
        end
    end
end

% 缓存根据不同A矩阵计算后的K，每行分别存储着时延为0,1,2,...,2d情况下的K
K_gather=zeros(2*d+1,P);
for i=1:(2*d+1)
    A=A_gather(:,(1:M)+(i-1)*M);
    K_gather(i,:)=[1,zeros(1,M-1)]*inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
end

% ********************************************************************************** %
% 第0组
% 不存在时延和丢包情况下的动态矩阵控制算法
for k=1:totalStep
    % 初始时刻，对相应值进行初始化
    if k==1
        % 控制量初始化
        U_control=zeros(1,simulationTime);
        U_control(:,1)=0;
        % 系统状态值初始化
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % 系统输出值初始化
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(:,1)=C_model_discrete*X_state(:,1);
        % 控制增量初始化
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % 初始预测值初始化
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
        % 实际使用控制增量初始化
        actualControlIncrement=zeros(1,simulationTime);
        % 校正向量初始化
        h=ones(N,1);
        % 输出误差初始化
        error=zeros(1,simulationTime);
        % Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        Y_setValue=ones(P,1)*Sv;
        continue
    end
    
    % 根据获得输出更新预测输出（校正），并移位
    error(k)=Y_outputValue(k)-Y_predictedValue(1,k);
    tempCorrection=Y_predictedValue(:,k)+correction_h*h*error(k);
    S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % 给出移位矩阵S
    Y_predictedValue(:,k)=S*tempCorrection;
    
    % 根据获得输出计算控制增量
    controlIncrement(:,k)=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,k));
    
    % 末尾时刻，跳出最后时刻的输出值计算，终止循环
    if k==totalStep
       continue; 
    end
    
    % 根据初始控制增量，计算下一时刻输出的初始预测值
    Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k);
    
    actualControlIncrement(k)=controlIncrement(1,k);
    % 发送控制增量后，计算实际控制量
    U_control(k)=U_control(k-1)+actualControlIncrement(k);
    
    X_state(:,k+1) = A_model_discrete * X_state(:,k) + B_model_discrete * U_control(k);
    Y_outputValue(k+1) = C_model_discrete * X_state(:,k+1);
    
end

% 实际使用控制量序列
figure(2)
stairs((1:totalStep)*T,U_control(1:totalStep),'-k');
grid on
hold on

% 系统输出曲线
figure(3)
plot((1:totalStep)*T,Y_outputValue(1:totalStep),'-k');
grid on
hold on

%{
% ********************************************************************************** %
% 第 1 组
% 问题：时延和乱序
% 方式：执行器，处理乱序（到达的数据包中选取时间戳最新的数据包）；控制器，处理乱序（主动丢包）
for k=1:totalStep
    % 初始时刻，对相应值进行初始化
    if k==1
        % 控制量初始化
        U_control=zeros(1,simulationTime);
        U_control(:,1)=0;
        % 系统状态值初始化
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % 系统输出值初始化
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(:,1)=C_model_discrete*X_state(:,1);
        % 控制增量初始化
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % 初始预测值初始化
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
        % 实际使用控制增量初始化
        actualControlValue=zeros(1,simulationTime);
        % 校正向量初始化
        h=ones(N,1);
        % 输出误差初始化
        error=zeros(1,simulationTime);
        % Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        Y_setValue=ones(P,1)*Sv;
        continue
    end
    
    % 根据获得输出更新预测输出（校正），并移位，考虑时延，控制器使用收到的最新输出值进行校正
    % 如果数据包乱序（主动丢包），跳过此次运算
    if timeSequenceInt(k)==(d+1)
        % 跳过乱序的数据包
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k);
        U_control(:,k)=U_control(:,k-1);
    else
        error(k)=Y_outputValue(k)-Y_predictedValue(1,k);
        tempCorrection=Y_predictedValue(:,k)+correction_h*h*error(k);
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % 给出移位矩阵S
        Y_predictedValue(:,k)=S*tempCorrection;
        % 根据获得输出计算控制增量，根据传感器-控制器的时延，选取合适的A矩阵计算控制增量
        % controlIncrement(:,k)=K_gather(1+timeSequenceInt_S2C(k-sigma_S2C(k)),:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        % 根据获得输出计算控制增量，不补偿时延
        controlIncrement(:,k)=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        % 根据初始控制增量，计算下一时刻输出的初始预测值
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k);
        
        U_control(:,k)=U_control(:,k-1)+controlIncrement(1,k);
    end
    
    % 末尾时刻，跳出最后时刻的输出值计算，终止循环
    if k==totalStep
       continue; 
    end
        
    % 发送控制量后，计算实际控制量
    actualControlValue(k)=U_control(k-sigma_SC_reorder(k));
    
    X_state(:,k+1) = A_model_discrete * X_state(:,k) + B_model_discrete * U_control(k);
    Y_outputValue(k+1) = C_model_discrete * X_state(:,k+1);
    
end

% 实际使用控制量序列
figure(2)
stairs((1:totalStep)*T,U_control(1:totalStep),'-.k');
grid on
hold on

% 系统输出曲线
figure(3)
plot((1:totalStep)*T,Y_outputValue(1:totalStep),'-.k');
grid on
hold on
%}

% ********************************************************************************** %
% 第 1 组
% 问题：时延和乱序
% 方式：执行器，处理乱序（到达的数据包中选取时间戳最新的数据包）；控制器，处理乱序（主动丢包）
for k=1:totalStep
    % 初始时刻，对相应值进行初始化
    if k==1
        % 控制量初始化
        U_control=zeros(1,simulationTime);
        U_control(:,1)=0;
        % 系统状态值初始化
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % 系统输出值初始化
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(:,1)=C_model_discrete*X_state(:,1);
        % 控制增量初始化
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % 初始预测值初始化
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
        % 实际使用控制增量初始化
        actualControlValue=zeros(1,simulationTime);
        % 校正向量初始化
        h=ones(N,1);
        % 输出误差初始化
        error=zeros(1,simulationTime);
        % Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        Y_setValue=ones(P,1)*Sv;
        continue
    end
    
    % 根据获得输出更新预测输出（校正），并移位，考虑时延，控制器使用收到的最新输出值进行校正
    % 如果数据包乱序（主动丢包），跳过此次运算
    if timeSequenceInt(k)==(d+1)
        % 跳过乱序的数据包
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k);
        U_control(:,k)=U_control(:,k-1);
    else
        error(k)=Y_outputValue(k)-Y_predictedValue(1,k);
        tempCorrection=Y_predictedValue(:,k)+correction_h*h*error(k);
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % 给出移位矩阵S
        Y_predictedValue(:,k)=S*tempCorrection;
        % 根据获得输出计算控制增量，根据传感器-控制器的时延，选取合适的A矩阵计算控制增量
        controlIncrement(:,k)=K_gather(1+timeSequenceInt(k),:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        % 根据获得输出计算控制增量，不补偿时延
        % controlIncrement(:,k)=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        % 根据初始控制增量，计算下一时刻输出的初始预测值
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k);
        
        U_control(:,k)=U_control(:,k-1)+controlIncrement(1,k);
    end
    
    % 末尾时刻，跳出最后时刻的输出值计算，终止循环
    if k==totalStep
       continue; 
    end
        
    % 发送控制量后，计算实际控制量
    actualControlValue(k)=U_control(k-sigma_SC_reorder(k));
    
    X_state(:,k+1) = A_model_discrete * X_state(:,k) + B_model_discrete * actualControlValue(k);
    Y_outputValue(k+1) = C_model_discrete * X_state(:,k+1);
    
end

% 实际使用控制量序列
figure(2)
stairs((1:totalStep)*T,U_control(1:totalStep),'-.k');
grid on
hold on

% 系统输出曲线
figure(3)
plot((1:totalStep)*T,Y_outputValue(1:totalStep),'-.k');
grid on
hold on

% ********************************************************************************** %
% 第 2 组
% 问题：时延和乱序
% 方式：执行器，不处理乱序；控制器，处理乱序（主动丢包）
% ********************************************************************************** %
% 第 3 组
% 问题：时延和乱序
% 方式：执行器，处理乱序（到达的数据包中选取时间戳最新的数据包）；控制器，不处理乱序（无主动丢包）
% ********************************************************************************** %
% 第 4 组
% 问题：时延和乱序
% 方式：执行器，不处理乱序（选取最新到达的数据包）；控制器，不处理乱序（无主动丢包）
for k=1:totalStep
    % 初始时刻，对相应值进行初始化
    if k==1
        index=k+1;
        % 控制量初始化
        U_control=zeros(1,simulationTime);
        U_control(:,1)=0;
        % 系统状态值初始化
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % 系统输出值初始化
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(:,1)=C_model_discrete*X_state(:,1);
        % 控制增量初始化
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % 初始预测值初始化
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
        % 实际使用控制增量初始化
        actualControlValue=zeros(1,simulationTime);
        % 校正向量初始化
        h=ones(N,1);
        % 输出误差初始化
        error=zeros(1,simulationTime);
        % Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        Y_setValue=ones(P,1)*Sv;
        continue
    end
    
    % 有数据包到达就开始控制量计算
    while index<=totalStep && timeStampList_S2C(index)<=k
        error(index)=Y_outputValue(timeStampList_S2C(index))-Y_predictedValue(1,index);
        tempCorrection=Y_predictedValue(:,index)+correction_h*h*error(index);
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % 给出移位矩阵S
        Y_predictedValue(:,index)=S*tempCorrection;

        % 根据获得输出计算控制增量，A矩阵不根据时延动态变化
        controlIncrement(:,timeStampList_S2C(index))=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,index));

        U_control(:,index)=U_control(:,index-1)+controlIncrement(1,index);
        % 根据初始控制增量，计算下一时刻输出的初始预测值
        Y_predictedValue(:,index+1)=Y_predictedValue(:,index)+stepResponse(:,1)*controlIncrement(1,index);
        
        index=index+1;
    end
    
    % 末尾时刻，跳出最后时刻的输出值计算，终止循环
    if k==totalStep
       continue; 
    end
    
    % 发送控制量后，计算控制器实际使用的控制量
    actualControlValue(k)=U_control(k-latestArrive_SC_noreorder(k));
    
    X_state(:,k+1) = A_model_discrete * X_state(:,k) + B_model_discrete * actualControlValue(k);
    Y_outputValue(k+1) = C_model_discrete * X_state(:,k+1);
    
end

% 实际使用控制量序列
figure(2)
set(gcf,'Position',[1000,400,700,500]);
stairs((1:totalStep)*T,U_control(1:totalStep),'--k');
legend('标准情况','处理乱序','不处理乱序');
xlabel('时间k');
ylabel('实际使用控制量u');
grid on
hold on

% 系统输出曲线
figure(3)
set(gcf,'Position',[1000,400,700,500]);
plot((1:totalStep)*T,Y_outputValue(1:totalStep),'--k');
legend('标准情况','处理乱序','不处理乱序');
xlabel('时间k');
ylabel('输出y');
grid on
hold on

