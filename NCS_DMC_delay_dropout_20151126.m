%  2015-11-26：参照NCS_DMC_delay_20151120.m中的程序，修改相对应的考虑丢包之后的情况；
%              用离散状态空间方程替代之前的纯粹基于阶跃响应系数的描述；
%              因为考虑了丢包，所以需要扩展A_gather矩阵，计算更多组(1+d+s)阶跃响应系数，s为最大连续丢包数；
%  2016-08-07：考在传感器和控制器之间存在乱序问题，将这一环节引入到代码中；
%  2016-08-12: 不考虑传感器和控制器之间的乱序问题，将问题简化成单边时延的情况；
% ********************************************************************************** %

clear all;clc;close all;

d=4; % 最大时延为 d*T
T=0.04; % 采样周期
totalStep=62;% 仿真总步长
% timeSequence= d*T*[zeros(1,50)]; % 无时延的情况
timeSequence = d*T*[0 rand(1,totalStep-3)-0.2 0 0];
timeSequenceInt = timeSequence;
timeSequenceLength = length(timeSequence);
controlSequence=zeros(2*(d+1),timeSequenceLength);

% 将时间序列转义为0,1,2,3......
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

% 设定丢包的概率，将发生丢包的时刻存入dropoutPoint，并将丢包时刻的时延设置为(d+1)T
dropoutRate=[0 rand(1,timeSequenceLength-3) 0 0];
dropoutPoint=[];
for i=1:timeSequenceLength
    if dropoutRate(i)>=0.9
        timeSequenceInt(i)=d+1;
        dropoutPoint=[dropoutPoint,i];
    end
end

% 排列成一个2(d+1)*timeSequenceLength的矩阵，用 第1,3,5,...,(2d+1)行，第k列 表示k时刻是否收到 k,(k-1),(k-2)...,(k-d)时刻的控制量
for i=1:timeSequenceLength
    for j=0:d
        if timeSequenceInt(i)==j
            controlSequence(2*j+1,i+j)=1;
        end
    end
end


disorderingPoint=[]; % 乱序点初始化
for k=1:timeSequenceLength
    
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
delayPoint=[];
for i=1:timeSequenceLength
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

%计算每个时刻的连续丢包数
continuousDropoutNum=zeros(1,timeSequenceLength);
for i=2:timeSequenceLength
    if i<=d+1
        dvalue=i-1;
    end
    for j=1:dvalue
        if timeSequenceInt(i-j)==(d+1)
            continuousDropoutNum(i)=continuousDropoutNum(i)+1;
        else
            break
        end
    end
end

figure(1)
t=T*(1:timeSequenceLength);
subplot(2,1,1);
% plot(t,timeSequenceInt*T,'o',t,timeSequence,'.');
plot(delayPoint*T,timeSequenceIntDraw*T,'*k',disorderingPoint*T,zeros(length(disorderingPoint),1),'ok',dropoutPoint*T,d*T,'xk');
% title('各个点的时延分布');
axis([0 timeSequenceLength*T -0.5*T (d+1)*T]);
legend('时延','发生乱序的时刻','丢包的数据');
xlabel('时间k');
ylabel('\tau _k');
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
M=10; % 控制时域
control_R=0.001; % 控制权矩阵系数
error_Q=1; % 误差权矩阵系数
correction_h=1; % 校正系数
Sv=1; % 设定值，参考轨迹
simulationTime=0+timeSequenceLength;

% 系统的离散状态空间方程
[A_model_discrete,B_model_discrete]=c2d(A_model,B_model,T);
C_model_discrete=C_model;

% 无时延系统的阶跃响应，获得阶跃响应系数
temp_response=dstep(A_model_discrete,B_model_discrete,C_model_discrete,0);
stepResponse(:,1)=[zeros(N,1),eye(N),zeros(N,length(temp_response)-(N+1))]*temp_response;

% 时延为delay=1:d的系统的阶跃响应系数
for delay=1:(2*d) 
    A_model_discrete_delay=expm(A_model*T);
    for i=1:N
        accumulation=zeros(2);%此处需要根据A_model、B_model的维数确定
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

% A_gather中依次排列着时延为0到dT情况下的动态矩阵
A_gather=zeros(P,M*(2*d+1));
for delay=0:(2*d)
    % 将阶跃响应系数赋给行向量a
    a=zeros(N,1);
    for i=1:N
        a(i)=stepResponse(i,delay+1); % a为系统阶跃响应的系数
    end
    % 计算不同时延情况下的动态矩阵
    for i=1:P
        for j=1:M
            if i-j+1>0
                A_gather(i,delay*M+j)=a(i-j+1);
            end
        end
    end
end

% 缓存根据不同A矩阵计算后的
K_gather=zeros(2*d+1,P);
for i=1:(2*d+1)
    A=A_gather(:,(1:M)+(i-1)*M);
    K_gather(i,:)=[1,zeros(1,M-1)]*inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
end

A=A_gather(:,1:M);



% ********************************************************************************** %
% 仅仅选取k+1时刻的最新控制量，若该时刻无控制量到达，保持上一时刻控制量，无reordering
% 反馈值，选择k+1时刻的最新输出值，若该时刻无输出值到达，保持上一时刻的输出值，无reordering
for k=1:timeSequenceLength-1 % N应该与timeSequenceLength保持一致
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        Y_setValue=ones(P,1)*Sv; 
        % 取最优控制量计算公式中可以离线计算的部分，将计算结果赋给K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % 控制量初始化
        controlValue=zeros(1,simulationTime);
        controlValue(1)=1;
        % 实际控制量初始化
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1);
        % 控制增量初始化
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % 校正向量
        h=ones(N,1); 
        % 输出误差初始化
        error=zeros(1,simulationTime);
        % 系统状态初始化
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % 输出值初始化
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % 实际获得输出值初始化
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % 输出预测值初始化
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
    
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k)); % 计算控制增量
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1); % 根据控制增量计算控制量
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % 根据控制增量计算(k+1)时刻y的预测值 

        if k<=d
            dValue=k;
        end
        % 根据时延，将控制量存储到controlSequence相应的位置
        for temp1=0:dValue
           if timeSequenceInt(k+1)==temp1
               controlSequence(2*(temp1+1),k+1+temp1)=controlValue(k+1);
           end
        end
        % 选取实际控制量actualControlValue
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
        % 施加控制量之后系统的实际输出和状态量
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % 根据(k+1)时刻y的实际值和预测值计算输出误差，以用于校正
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % 运用启发式误差预测方法，修正(k+1)时刻的预测值
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % 给出移位矩阵S
        Y_predictedValue(:,k+1)=S*tempCorrection; % 给出(k+1)时刻
    
end

% ********************************************************************************** %
% 实际使用控制量序列
figure(2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'-.k');
grid on
hold on

% 系统输出曲线
figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-.k');
grid on
hold on

% ********************************************************************************** %
% 处理了产生disordering的控制量，有reordering
for k=1:timeSequenceLength-1 % N应该与timeSequenceLength保持一致
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        Y_setValue=ones(P,1)*Sv; 
        % 取最优控制量计算公式中可以离线计算的部分，将计算结果赋给K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % 控制量初始化
        controlValue=zeros(1,simulationTime);
        controlValue(1)=1;
        % 实际控制量初始化
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1);
        % 控制增量初始化
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % 校正向量
        h=ones(N,1); 
        % 输出误差初始化
        error=zeros(1,simulationTime);
        % 系统状态初始化
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % 输出值初始化
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % 实际获得输出值初始化
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % 输出预测值初始化
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
    
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k)); % 计算控制增量
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1); % 根据控制增量计算控制量
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % 根据控制增量计算(k+1)时刻y的预测值 

        if k<=d+1
            dValue=k;
        end
        temp3=zeros(dValue+1,1);
        for temp2=0:dValue
            temp3(temp2+1)=timeSequenceInt(k+1-temp2)-(temp2+continuousDropoutNum(k+1-temp2));
        end
        % minIndex,最小值索引
        % minValue，最小值
        for temp2=0:dValue
            if temp3(temp2+1)<=0
                minIndex=temp2+continuousDropoutNum(k+1-temp2);
                minValue=temp3(minIndex+1);
                break
            end
        end
        
        sigma(k+1)=minIndex;
        actualControlValue(k+1)=controlValue(k+1-minIndex);
%         actualControlValue(k+1)=actualControlValue(k)+controlIncrement(1,k+1-minIndex);

        % 施加控制量之后系统的实际输出和状态量
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % 根据(k+1)时刻y的实际值和预测值计算误差，以用于校正
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % 运用启发式误差预测方法，修正(k+1)时刻的预测值
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % 给出移位矩阵S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

figure(1)
subplot(2,1,2);
plot((1:timeSequenceLength)*T,sigma,'*k')
xlabel('时间k');
ylabel('\sigma(k)');
grid on
hold on

% 实际使用控制量序列
figure(2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'--k');
% title('控制量');
xlabel('Time/s');
ylabel('控制量输出');
grid on
hold on

% 系统输出曲线
figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'--k');
% title('输出曲线');
xlabel('时间k');
ylabel('输出y');
grid on
hold on

% ********************************************************************************** %
% 根据时延，准确补偿时延和丢包对系统造成的影响
for k=1:timeSequenceLength-1 % N应该与timeSequenceLength保持一致
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        Y_setValue=ones(P,1)*Sv; 
        % 取最优控制量计算公式中可以离线计算的部分，将计算结果赋给K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % 控制量初始化
        controlValue=zeros(1,simulationTime);
        controlValue(1)=1;
        % 实际控制量初始化
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1);
        % 控制增量初始化
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % 校正向量
        h=ones(N,1); 
        % 输出误差初始化
        error=zeros(1,simulationTime);
        % 系统状态初始化
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % 输出值初始化
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % 实际获得输出值初始化
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % 输出预测值初始化
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
        
        A=A_gather(:,((M*timeSequenceInt(k+1)+1):(M*timeSequenceInt(k+1)+M)));
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k)); % 计算控制增量
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1); % 根据控制增量计算控制量
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % 根据控制增量计算(k+1)时刻y的预测值 
        
        minIndex=sigma(k+1);
        actualControlValue(k+1)=controlValue(k+1-minIndex);
%         actualControlValue(k+1)=actualControlValue(k)+controlIncrement(1,k+1-minIndex);

        % 施加控制量之后系统的实际输出和状态量
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % 根据(k+1)时刻y的实际值和预测值计算误差，以用于校正
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % 运用启发式误差预测方法，修正(k+1)时刻的预测值
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % 给出移位矩阵S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

% ********************************************************************************** %
% 实际使用控制量序列
figure(2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'-k');
legend('不处理乱序问题','处理了乱序问题','改进DMC后');
xlabel('Time/s');
ylabel('控制量输出');
grid on
hold on

% 系统输出曲线
figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-k');
legend('不处理乱序问题','处理了乱序问题','改进DMC后');
xlabel('时间k');
ylabel('输出y');
grid on
hold on

%{
% ********************************************************************************** %
% 根据时延，准确补偿时延和丢包对系统造成的影响
for k=1:timeSequenceLength-1 % N应该与timeSequenceLength保持一致
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        Y_setValue=ones(P,1)*Sv; 
        % 取最优控制量计算公式中可以离线计算的部分，将计算结果赋给K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % 控制量初始化
        controlValue=zeros(1+2*d,simulationTime);
        controlValue(1,1)=1;
        % 实际控制量初始化
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1,1);
        % 控制增量初始化
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % 校正向量
        h=ones(N,1); 
        % 输出误差初始化
        error=zeros(1,simulationTime);
        % 系统状态初始化
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % 输出值初始化
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % 实际获得输出值初始化
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % 输出预测值初始化
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
        
%         A=A_gather(:,((M*timeSequenceInt(k+1)+1):(M*timeSequenceInt(k+1)+M)));
%         K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        
        controlIncrement(1,k+1)=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        for i=1:(2*d+1)
            controlValue(i,k+1)=controlValue(1,k)+K_gather(i,:)*(Y_setValue-Y_predictedValue(1:P,k));
        end
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % 根据控制增量计算(k+1)时刻y的预测值 
        
        if k<=d+1
            dValue=k;
        end
        temp3=zeros(dValue+1,1);
        for temp2=0:dValue
            temp3(temp2+1)=timeSequenceInt(k+1-temp2)-(temp2+continuousDropoutNum(k+1-temp2));
        end
        
        minIndex=sigma(k+1);
        actualControlValue(k+1)=controlValue(1,k+1-minIndex);
%         actualControlValue(k+1)=actualControlValue(k)+controlIncrement(1,k+1-minIndex);

        % 施加控制量之后系统的实际输出和状态量
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % 根据(k+1)时刻y的实际值和预测值计算误差，以用于校正
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % 运用启发式误差预测方法，修正(k+1)时刻的预测值
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % 给出移位矩阵S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-');
hold on
%}

% ********************************************************************************** %
% 使用冗余的控制增量来代替丢包或延时的控制量
%{
for k=1:timeSequenceLength-1 % N应该与timeSequenceLength保持一致
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',参考轨迹，即期望值
        Y_setValue=ones(P,1)*Sv; 
        % 取最优控制量计算公式中可以离线计算的部分，将计算结果赋给K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % 控制量初始化
        controlValue=zeros(1,simulationTime);
        controlValue(1)=1;
        % 实际控制量初始化
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1);
        % 控制增量初始化
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % 校正向量
        h=ones(N,1); 
        % 输出误差初始化
        error=zeros(1,simulationTime);
        % 系统状态初始化
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % 输出值初始化
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % 实际获得输出值初始化
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % 输出预测值初始化
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
        
        A=A_gather(:,1:M);
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k)); % 计算控制增量
        controlValue(k+1)=controlValue(k)+controlIncrement(1+timeSequenceInt(k+1),k+1); % 根据控制增量计算控制量
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % 根据控制增量计算(k+1)时刻y的预测值 
        
        if k<=d+1
            dValue=k;
        end
        temp3=zeros(dValue+1,1);
        for temp2=0:dValue
            temp3(temp2+1)=timeSequenceInt(k+1-temp2)-(temp2+continuousDropoutNum(k+1-temp2));
        end
        % minIndex,最小值索引
        % minValue，最小值
        for temp2=0:dValue
            if temp3(temp2+1)<=0
                minIndex=temp2+continuousDropoutNum(k+1-temp2);
                minValue=temp3(minIndex+1);
                break
            end
        end
        
        sigma(k+1)=minIndex;
        actualControlValue(k+1)=controlValue(k+1-minIndex);
%         actualControlValue(k+1)=actualControlValue(k)+controlIncrement(1,k+1-minIndex);

        % 施加控制量之后系统的实际输出和状态量
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % 根据(k+1)时刻y的实际值和预测值计算误差，以用于校正
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % 运用启发式误差预测方法，修正(k+1)时刻的预测值
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % 给出移位矩阵S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end
figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-');
hold on
%}