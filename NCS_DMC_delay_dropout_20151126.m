%  2015-11-26������NCS_DMC_delay_20151120.m�еĳ����޸����Ӧ�Ŀ��Ƕ���֮��������
%              ����ɢ״̬�ռ䷽�����֮ǰ�Ĵ�����ڽ�Ծ��Ӧϵ����������
%              ��Ϊ�����˶�����������Ҫ��չA_gather���󣬼��������(1+d+s)��Ծ��Ӧϵ����sΪ���������������
%  2016-08-07�����ڴ������Ϳ�����֮������������⣬����һ�������뵽�����У�
%  2016-08-12: �����Ǵ������Ϳ�����֮����������⣬������򻯳ɵ���ʱ�ӵ������
% ********************************************************************************** %

clear all;clc;close all;

d=4; % ���ʱ��Ϊ d*T
T=0.04; % ��������
totalStep=62;% �����ܲ���
% timeSequence= d*T*[zeros(1,50)]; % ��ʱ�ӵ����
timeSequence = d*T*[0 rand(1,totalStep-3)-0.2 0 0];
timeSequenceInt = timeSequence;
timeSequenceLength = length(timeSequence);
controlSequence=zeros(2*(d+1),timeSequenceLength);

% ��ʱ������ת��Ϊ0,1,2,3......
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

% �趨�����ĸ��ʣ�������������ʱ�̴���dropoutPoint����������ʱ�̵�ʱ������Ϊ(d+1)T
dropoutRate=[0 rand(1,timeSequenceLength-3) 0 0];
dropoutPoint=[];
for i=1:timeSequenceLength
    if dropoutRate(i)>=0.9
        timeSequenceInt(i)=d+1;
        dropoutPoint=[dropoutPoint,i];
    end
end

% ���г�һ��2(d+1)*timeSequenceLength�ľ����� ��1,3,5,...,(2d+1)�У���k�� ��ʾkʱ���Ƿ��յ� k,(k-1),(k-2)...,(k-d)ʱ�̵Ŀ�����
for i=1:timeSequenceLength
    for j=0:d
        if timeSequenceInt(i)==j
            controlSequence(2*j+1,i+j)=1;
        end
    end
end


disorderingPoint=[]; % ������ʼ��
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

%ȡ��δ����������ʱ�̵�ʱ�ӣ�����timeSequenceIntDraw�У�����δ������ʱ�̴���delayPoint��
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

%����ÿ��ʱ�̵�����������
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
% title('�������ʱ�ӷֲ�');
axis([0 timeSequenceLength*T -0.5*T (d+1)*T]);
legend('ʱ��','���������ʱ��','����������');
xlabel('ʱ��k');
ylabel('\tau _k');
grid on
hold on

% ********************************************************************************** %
% ���ض����״̬�ռ䷽��
A_model=[-4 -0.03;0.75 -10];
B_model=[2 0]';
C_model=[0 1];

% Ԥ����DMC����������Ҫ����ز���
N=40; % ģ��ʱ��
P=15; % Ԥ��ʱ��
M=10; % ����ʱ��
control_R=0.001; % ����Ȩ����ϵ��
error_Q=1; % ���Ȩ����ϵ��
correction_h=1; % У��ϵ��
Sv=1; % �趨ֵ���ο��켣
simulationTime=0+timeSequenceLength;

% ϵͳ����ɢ״̬�ռ䷽��
[A_model_discrete,B_model_discrete]=c2d(A_model,B_model,T);
C_model_discrete=C_model;

% ��ʱ��ϵͳ�Ľ�Ծ��Ӧ����ý�Ծ��Ӧϵ��
temp_response=dstep(A_model_discrete,B_model_discrete,C_model_discrete,0);
stepResponse(:,1)=[zeros(N,1),eye(N),zeros(N,length(temp_response)-(N+1))]*temp_response;

% ʱ��Ϊdelay=1:d��ϵͳ�Ľ�Ծ��Ӧϵ��
for delay=1:(2*d) 
    A_model_discrete_delay=expm(A_model*T);
    for i=1:N
        accumulation=zeros(2);%�˴���Ҫ����A_model��B_model��ά��ȷ��
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

% A_gather������������ʱ��Ϊ0��dT����µĶ�̬����
A_gather=zeros(P,M*(2*d+1));
for delay=0:(2*d)
    % ����Ծ��Ӧϵ������������a
    a=zeros(N,1);
    for i=1:N
        a(i)=stepResponse(i,delay+1); % aΪϵͳ��Ծ��Ӧ��ϵ��
    end
    % ���㲻ͬʱ������µĶ�̬����
    for i=1:P
        for j=1:M
            if i-j+1>0
                A_gather(i,delay*M+j)=a(i-j+1);
            end
        end
    end
end

% ������ݲ�ͬA���������
K_gather=zeros(2*d+1,P);
for i=1:(2*d+1)
    A=A_gather(:,(1:M)+(i-1)*M);
    K_gather(i,:)=[1,zeros(1,M-1)]*inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
end

A=A_gather(:,1:M);



% ********************************************************************************** %
% ����ѡȡk+1ʱ�̵����¿�����������ʱ���޿��������������һʱ�̿���������reordering
% ����ֵ��ѡ��k+1ʱ�̵��������ֵ������ʱ�������ֵ���������һʱ�̵����ֵ����reordering
for k=1:timeSequenceLength-1 % NӦ����timeSequenceLength����һ��
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        Y_setValue=ones(P,1)*Sv; 
        % ȡ���ſ��������㹫ʽ�п������߼���Ĳ��֣�������������K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % ��������ʼ��
        controlValue=zeros(1,simulationTime);
        controlValue(1)=1;
        % ʵ�ʿ�������ʼ��
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1);
        % ����������ʼ��
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % У������
        h=ones(N,1); 
        % �������ʼ��
        error=zeros(1,simulationTime);
        % ϵͳ״̬��ʼ��
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % ���ֵ��ʼ��
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % ʵ�ʻ�����ֵ��ʼ��
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % ���Ԥ��ֵ��ʼ��
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
    
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k)); % �����������
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1); % ���ݿ����������������
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % ���ݿ�����������(k+1)ʱ��y��Ԥ��ֵ 

        if k<=d
            dValue=k;
        end
        % ����ʱ�ӣ����������洢��controlSequence��Ӧ��λ��
        for temp1=0:dValue
           if timeSequenceInt(k+1)==temp1
               controlSequence(2*(temp1+1),k+1+temp1)=controlValue(k+1);
           end
        end
        % ѡȡʵ�ʿ�����actualControlValue
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
        % ʩ�ӿ�����֮��ϵͳ��ʵ�������״̬��
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % ����(k+1)ʱ��y��ʵ��ֵ��Ԥ��ֵ���������������У��
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % ��������ʽ���Ԥ�ⷽ��������(k+1)ʱ�̵�Ԥ��ֵ
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % ������λ����S
        Y_predictedValue(:,k+1)=S*tempCorrection; % ����(k+1)ʱ��
    
end

% ********************************************************************************** %
% ʵ��ʹ�ÿ���������
figure(2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'-.k');
grid on
hold on

% ϵͳ�������
figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-.k');
grid on
hold on

% ********************************************************************************** %
% �����˲���disordering�Ŀ���������reordering
for k=1:timeSequenceLength-1 % NӦ����timeSequenceLength����һ��
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        Y_setValue=ones(P,1)*Sv; 
        % ȡ���ſ��������㹫ʽ�п������߼���Ĳ��֣�������������K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % ��������ʼ��
        controlValue=zeros(1,simulationTime);
        controlValue(1)=1;
        % ʵ�ʿ�������ʼ��
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1);
        % ����������ʼ��
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % У������
        h=ones(N,1); 
        % �������ʼ��
        error=zeros(1,simulationTime);
        % ϵͳ״̬��ʼ��
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % ���ֵ��ʼ��
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % ʵ�ʻ�����ֵ��ʼ��
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % ���Ԥ��ֵ��ʼ��
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
    
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k)); % �����������
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1); % ���ݿ����������������
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % ���ݿ�����������(k+1)ʱ��y��Ԥ��ֵ 

        if k<=d+1
            dValue=k;
        end
        temp3=zeros(dValue+1,1);
        for temp2=0:dValue
            temp3(temp2+1)=timeSequenceInt(k+1-temp2)-(temp2+continuousDropoutNum(k+1-temp2));
        end
        % minIndex,��Сֵ����
        % minValue����Сֵ
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

        % ʩ�ӿ�����֮��ϵͳ��ʵ�������״̬��
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % ����(k+1)ʱ��y��ʵ��ֵ��Ԥ��ֵ������������У��
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % ��������ʽ���Ԥ�ⷽ��������(k+1)ʱ�̵�Ԥ��ֵ
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % ������λ����S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

figure(1)
subplot(2,1,2);
plot((1:timeSequenceLength)*T,sigma,'*k')
xlabel('ʱ��k');
ylabel('\sigma(k)');
grid on
hold on

% ʵ��ʹ�ÿ���������
figure(2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'--k');
% title('������');
xlabel('Time/s');
ylabel('���������');
grid on
hold on

% ϵͳ�������
figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'--k');
% title('�������');
xlabel('ʱ��k');
ylabel('���y');
grid on
hold on

% ********************************************************************************** %
% ����ʱ�ӣ�׼ȷ����ʱ�ӺͶ�����ϵͳ��ɵ�Ӱ��
for k=1:timeSequenceLength-1 % NӦ����timeSequenceLength����һ��
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        Y_setValue=ones(P,1)*Sv; 
        % ȡ���ſ��������㹫ʽ�п������߼���Ĳ��֣�������������K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % ��������ʼ��
        controlValue=zeros(1,simulationTime);
        controlValue(1)=1;
        % ʵ�ʿ�������ʼ��
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1);
        % ����������ʼ��
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % У������
        h=ones(N,1); 
        % �������ʼ��
        error=zeros(1,simulationTime);
        % ϵͳ״̬��ʼ��
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % ���ֵ��ʼ��
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % ʵ�ʻ�����ֵ��ʼ��
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % ���Ԥ��ֵ��ʼ��
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
        
        A=A_gather(:,((M*timeSequenceInt(k+1)+1):(M*timeSequenceInt(k+1)+M)));
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k)); % �����������
        controlValue(k+1)=controlValue(k)+controlIncrement(1,k+1); % ���ݿ����������������
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % ���ݿ�����������(k+1)ʱ��y��Ԥ��ֵ 
        
        minIndex=sigma(k+1);
        actualControlValue(k+1)=controlValue(k+1-minIndex);
%         actualControlValue(k+1)=actualControlValue(k)+controlIncrement(1,k+1-minIndex);

        % ʩ�ӿ�����֮��ϵͳ��ʵ�������״̬��
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % ����(k+1)ʱ��y��ʵ��ֵ��Ԥ��ֵ������������У��
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % ��������ʽ���Ԥ�ⷽ��������(k+1)ʱ�̵�Ԥ��ֵ
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % ������λ����S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

% ********************************************************************************** %
% ʵ��ʹ�ÿ���������
figure(2)
stairs((1:timeSequenceLength)*T,actualControlValue(1:timeSequenceLength),'-k');
legend('��������������','��������������','�Ľ�DMC��');
xlabel('Time/s');
ylabel('���������');
grid on
hold on

% ϵͳ�������
figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-k');
legend('��������������','��������������','�Ľ�DMC��');
xlabel('ʱ��k');
ylabel('���y');
grid on
hold on

%{
% ********************************************************************************** %
% ����ʱ�ӣ�׼ȷ����ʱ�ӺͶ�����ϵͳ��ɵ�Ӱ��
for k=1:timeSequenceLength-1 % NӦ����timeSequenceLength����һ��
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        Y_setValue=ones(P,1)*Sv; 
        % ȡ���ſ��������㹫ʽ�п������߼���Ĳ��֣�������������K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % ��������ʼ��
        controlValue=zeros(1+2*d,simulationTime);
        controlValue(1,1)=1;
        % ʵ�ʿ�������ʼ��
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1,1);
        % ����������ʼ��
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % У������
        h=ones(N,1); 
        % �������ʼ��
        error=zeros(1,simulationTime);
        % ϵͳ״̬��ʼ��
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % ���ֵ��ʼ��
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % ʵ�ʻ�����ֵ��ʼ��
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % ���Ԥ��ֵ��ʼ��
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
        
%         A=A_gather(:,((M*timeSequenceInt(k+1)+1):(M*timeSequenceInt(k+1)+M)));
%         K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        
        controlIncrement(1,k+1)=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        for i=1:(2*d+1)
            controlValue(i,k+1)=controlValue(1,k)+K_gather(i,:)*(Y_setValue-Y_predictedValue(1:P,k));
        end
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % ���ݿ�����������(k+1)ʱ��y��Ԥ��ֵ 
        
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

        % ʩ�ӿ�����֮��ϵͳ��ʵ�������״̬��
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % ����(k+1)ʱ��y��ʵ��ֵ��Ԥ��ֵ������������У��
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % ��������ʽ���Ԥ�ⷽ��������(k+1)ʱ�̵�Ԥ��ֵ
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % ������λ����S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end

figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-');
hold on
%}

% ********************************************************************************** %
% ʹ������Ŀ������������涪������ʱ�Ŀ�����
%{
for k=1:timeSequenceLength-1 % NӦ����timeSequenceLength����һ��
    
    if k==1
        % Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        Y_setValue=ones(P,1)*Sv; 
        % ȡ���ſ��������㹫ʽ�п������߼���Ĳ��֣�������������K
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        % ��������ʼ��
        controlValue=zeros(1,simulationTime);
        controlValue(1)=1;
        % ʵ�ʿ�������ʼ��
        actualControlValue=zeros(1,simulationTime);
        actualControlValue=controlValue(1);
        % ����������ʼ��
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % У������
        h=ones(N,1); 
        % �������ʼ��
        error=zeros(1,simulationTime);
        % ϵͳ״̬��ʼ��
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % ���ֵ��ʼ��
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(1)=C_model_discrete*X_state(:,1);
        % ʵ�ʻ�����ֵ��ʼ��
        actualOutputValue=zeros(1,simulationTime);
        actualOutputValue=Y_outputValue(1);
        % ���Ԥ��ֵ��ʼ��
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
    end
        
        A=A_gather(:,1:M);
        K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
        controlIncrement(:,k+1)=K*(Y_setValue-Y_predictedValue(1:P,k)); % �����������
        controlValue(k+1)=controlValue(k)+controlIncrement(1+timeSequenceInt(k+1),k+1); % ���ݿ����������������
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k+1); % ���ݿ�����������(k+1)ʱ��y��Ԥ��ֵ 
        
        if k<=d+1
            dValue=k;
        end
        temp3=zeros(dValue+1,1);
        for temp2=0:dValue
            temp3(temp2+1)=timeSequenceInt(k+1-temp2)-(temp2+continuousDropoutNum(k+1-temp2));
        end
        % minIndex,��Сֵ����
        % minValue����Сֵ
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

        % ʩ�ӿ�����֮��ϵͳ��ʵ�������״̬��
        X_state(:,k+1)=A_model_discrete*X_state(:,k)+B_model_discrete*actualControlValue(k+1);
        Y_outputValue(k+1)=C_model_discrete*X_state(:,k+1);
        
        error(k+1)=Y_outputValue(k+1)-Y_predictedValue(1,k+1); % ����(k+1)ʱ��y��ʵ��ֵ��Ԥ��ֵ������������У��
        tempCorrection=Y_predictedValue(:,k+1)+correction_h*h*error(k+1); % ��������ʽ���Ԥ�ⷽ��������(k+1)ʱ�̵�Ԥ��ֵ
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % ������λ����S
        Y_predictedValue(:,k+1)=S*tempCorrection;
    
end
figure(3)
plot((1:timeSequenceLength)*T,Y_outputValue(1:timeSequenceLength),'-');
hold on
%}