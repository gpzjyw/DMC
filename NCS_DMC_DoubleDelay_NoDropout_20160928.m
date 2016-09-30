%  2016-08-23�����������¸Ľ���1.����˫��ʱ�ӵ�������Ҷ����ܳ�����������⣻
%                             2.����Ĵ���ִ����ѡȡ���¿����źţ�������Ҳѡ�����µ����ֵ��
%                             3.��������У�����֣����õĿ�����������ʱ�� ������-������ ʱ��ȷ����
%                             4.�����������¼��������յ�һ�����ֵ����һ���¼�����ѡ�õ����ֵΪ���µ����ֵ��
%  2016-09-04�����������¿���: 1.ȥ���������أ���������ʱ�Ӻ�����������
%                             2.��������ִ�������û�����������ȷ����ʱ��ת��Ϊ���ɹ̶�ʱ�� 0,1,2,3,...,d��
%  2016-09-20�����������¿���: 1.�����������¼�������ִ��������ʱ�������Ĺ�����ʽ��
%                             2.����������������ݰ���������������
%  2016-09-28��������һ�����飺1.�����������⣬������ʱ�ӽ��в���
%                             2.����Ϊ���������
% ********************************************************************************** %

clear all;clc;close all;

d=4; % ���ʱ��Ϊ d*T
T=0.04; % ��������
totalStep=75;% �����ܲ���
timeSequence = d*T*[0 rand(1,totalStep-d-1)-0.2 zeros(1,d)]; % ÿ��ʱ�����ݰ���ʱ�Ӵ�С������ʱ�ӵĸ���Ϊ80%��
timeSequenceInt = timeSequence; % �������Ĵ��ڣ���ʱ��ת��Ϊ�̶���������ʱ�ӣ�0,T,2T,...,dT

% ��ʱ��
% ��ʱ������ת��Ϊ0,1,2,3....
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

% ����sigma(k)��ִ�������ݸ�ֵѡȡ���µĿ������������ڣ�����������������ִ������������
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

% ����ʱ�ӣ�������һʱ�����µ���Ŀ������������ڣ�����������������ִ��������������
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

% ������-������ʱ��
% �� ������-������ʱ�� С�ڵ��� ȫ��ʱ�� ΪԼ�������� ������-������ʱ��
timeSequence_S2C=zeros(1,totalStep);
timeSequenceInt_S2C=zeros(1,totalStep);
for i=1:totalStep
    timeSequence_S2C(i) = getRandom(0, timeSequence(i));
    timeSequenceInt_S2C(i) = ceil(timeSequence_S2C(i)/T);
end

% ȷ�� ������-������ ��������ݰ��������ʱ�����ݰ�Ϊ��������ʱ��d+1��ʾ��
disorderingPoint_S2C=[]; % ����ʱ�̵ļ���
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

% ���ݾ��嵽��ʱ�����������
arriveTime_S2C=[[1:totalStep] ; [1:totalStep]*T+timeSequence_S2C];
arriveTimeOrder_S2C=sortrows(arriveTime_S2C', 2)';
% ���������յ����ݰ�ʱ��˳��
timeStampList_S2C=arriveTimeOrder_S2C(1,:);
% ÿ������ʱ�̵������ݰ���ʱ������������°�˳������
timeStampSquare_S2C=zeros(d+1,totalStep);

% ����ʱ��ȷ���������ݰ��ļ���
disorderingPoint=[]; % ����ʱ�̵ļ���
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

% δ���������ݰ���ʱ�ӣ�����timeSequenceIntDraw
% δ���������ݰ���ʱ���������delayPoint
% ���������ݰ���ʱ���������dropoutPoint
timeSequenceIntDraw=[];
delayPoint=[]; % ʱ��ʱ�̵ļ���
dropoutPoint=[]; % ����ʱ�̼���
for i=1:totalStep
    if timeSequenceInt(i)==(d+1)
        dropoutPoint=[dropoutPoint,i];
        continue;
    end
    timeSequenceIntDraw=[timeSequenceIntDraw,timeSequenceInt(i)];
    delayPoint=[delayPoint,i];
end

% ����sigma(k)��ִ�������ݸ�ֵѡȡ���µĿ��������������ڣ���������������ִ������������
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

% ����ʱ�ӣ�������һʱ�����µ���Ŀ��������������ڣ���������������ִ��������������
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

% ��ͼ
figure(1)
set(gcf,'Position',[100,400,700,500]);
t=T*(1:totalStep);
subplot(2,1,1);
plot(delayPoint*T,timeSequenceIntDraw*T,'*k',disorderingPoint*T,zeros(length(disorderingPoint),1),'ok',dropoutPoint*T,d*T,'xk');
axis([0 totalStep*T -0.5*T (d+1)*T]);
legend('���ݰ���ʱ��','ִ�������������ʱ��','�������������ݰ�','Orientation','horizontal','Location','best');
xlabel('ʱ��k');
ylabel('\tau _k');
grid on
hold on
subplot(2,1,2);
plot((1:totalStep)*T,sigma_SC_reorder,'*k')
xlabel('ʱ��k');
ylabel('\sigma(k)');
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
M=5; % ����ʱ��
control_R=0.001; % ����Ȩ����ϵ��
error_Q=1; % ���Ȩ����ϵ��
correction_h=1; % У��ϵ��
Sv=1; % �趨ֵ���ο��켣
simulationTime=0+totalStep;

% ϵͳ����ɢ״̬�ռ䷽��
[A_model_discrete,B_model_discrete]=c2d(A_model,B_model,T);
C_model_discrete=C_model;

% ��ʱ��ϵͳ�Ľ�Ծ��Ӧ����ý�Ծ��Ӧϵ��
temp_response=dstep(A_model_discrete,B_model_discrete,C_model_discrete,0);
stepResponse(:,1)=[zeros(N,1),eye(N),zeros(N,length(temp_response)-(N+1))]*temp_response;

% ʱ��Ϊ1,2,3,...,dʱ��Ľ�Ծ��Ӧϵ��
% stepResponse ��ÿ�д洢��ʱ��Ϊ0,1,2,...,2d����µĽ�Ծ��Ӧϵ��
for delay=1:(2*d)
    for i=1:N
        accumulation=zeros(2,2);%�˴���Ҫ����A_model��B_model��ά��ȷ��
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

% A_gather������������ʱ��Ϊ0��2d�µĶ�̬����
A_gather=zeros(P,M*(2*d+1));
for delay=0:(2*d)
    % ����Ծ��Ӧϵ������������a
    a=zeros(N,1);
    a=stepResponse(:,delay+1); % aΪϵͳ��Ծ��Ӧ��ϵ��
    % ���㲻ͬʱ������µĶ�̬����
    for i=1:P
        for j=1:M
            if i-j+1>0
                A_gather(i,delay*M+j)=a(i-j+1);
            end
        end
    end
end

% ������ݲ�ͬA���������K��ÿ�зֱ�洢��ʱ��Ϊ0,1,2,...,2d����µ�K
K_gather=zeros(2*d+1,P);
for i=1:(2*d+1)
    A=A_gather(:,(1:M)+(i-1)*M);
    K_gather(i,:)=[1,zeros(1,M-1)]*inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A'*eye(P)*error_Q;
end

% ********************************************************************************** %
% ��0��
% ������ʱ�ӺͶ�������µĶ�̬��������㷨
for k=1:totalStep
    % ��ʼʱ�̣�����Ӧֵ���г�ʼ��
    if k==1
        % ��������ʼ��
        U_control=zeros(1,simulationTime);
        U_control(:,1)=0;
        % ϵͳ״ֵ̬��ʼ��
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % ϵͳ���ֵ��ʼ��
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(:,1)=C_model_discrete*X_state(:,1);
        % ����������ʼ��
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % ��ʼԤ��ֵ��ʼ��
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
        % ʵ��ʹ�ÿ���������ʼ��
        actualControlIncrement=zeros(1,simulationTime);
        % У��������ʼ��
        h=ones(N,1);
        % �������ʼ��
        error=zeros(1,simulationTime);
        % Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        Y_setValue=ones(P,1)*Sv;
        continue
    end
    
    % ���ݻ���������Ԥ�������У����������λ
    error(k)=Y_outputValue(k)-Y_predictedValue(1,k);
    tempCorrection=Y_predictedValue(:,k)+correction_h*h*error(k);
    S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % ������λ����S
    Y_predictedValue(:,k)=S*tempCorrection;
    
    % ���ݻ����������������
    controlIncrement(:,k)=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,k));
    
    % ĩβʱ�̣��������ʱ�̵����ֵ���㣬��ֹѭ��
    if k==totalStep
       continue; 
    end
    
    % ���ݳ�ʼ����������������һʱ������ĳ�ʼԤ��ֵ
    Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k);
    
    actualControlIncrement(k)=controlIncrement(1,k);
    % ���Ϳ��������󣬼���ʵ�ʿ�����
    U_control(k)=U_control(k-1)+actualControlIncrement(k);
    
    X_state(:,k+1) = A_model_discrete * X_state(:,k) + B_model_discrete * U_control(k);
    Y_outputValue(k+1) = C_model_discrete * X_state(:,k+1);
    
end

% ʵ��ʹ�ÿ���������
figure(2)
stairs((1:totalStep)*T,U_control(1:totalStep),'-k');
grid on
hold on

% ϵͳ�������
figure(3)
plot((1:totalStep)*T,Y_outputValue(1:totalStep),'-k');
grid on
hold on

%{
% ********************************************************************************** %
% �� 1 ��
% ���⣺ʱ�Ӻ�����
% ��ʽ��ִ�������������򣨵�������ݰ���ѡȡʱ������µ����ݰ�������������������������������
for k=1:totalStep
    % ��ʼʱ�̣�����Ӧֵ���г�ʼ��
    if k==1
        % ��������ʼ��
        U_control=zeros(1,simulationTime);
        U_control(:,1)=0;
        % ϵͳ״ֵ̬��ʼ��
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % ϵͳ���ֵ��ʼ��
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(:,1)=C_model_discrete*X_state(:,1);
        % ����������ʼ��
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % ��ʼԤ��ֵ��ʼ��
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
        % ʵ��ʹ�ÿ���������ʼ��
        actualControlValue=zeros(1,simulationTime);
        % У��������ʼ��
        h=ones(N,1);
        % �������ʼ��
        error=zeros(1,simulationTime);
        % Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        Y_setValue=ones(P,1)*Sv;
        continue
    end
    
    % ���ݻ���������Ԥ�������У����������λ������ʱ�ӣ�������ʹ���յ����������ֵ����У��
    % ������ݰ����������������������˴�����
    if timeSequenceInt(k)==(d+1)
        % ������������ݰ�
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k);
        U_control(:,k)=U_control(:,k-1);
    else
        error(k)=Y_outputValue(k)-Y_predictedValue(1,k);
        tempCorrection=Y_predictedValue(:,k)+correction_h*h*error(k);
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % ������λ����S
        Y_predictedValue(:,k)=S*tempCorrection;
        % ���ݻ���������������������ݴ�����-��������ʱ�ӣ�ѡȡ���ʵ�A��������������
        % controlIncrement(:,k)=K_gather(1+timeSequenceInt_S2C(k-sigma_S2C(k)),:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        % ���ݻ������������������������ʱ��
        controlIncrement(:,k)=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        % ���ݳ�ʼ����������������һʱ������ĳ�ʼԤ��ֵ
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k);
        
        U_control(:,k)=U_control(:,k-1)+controlIncrement(1,k);
    end
    
    % ĩβʱ�̣��������ʱ�̵����ֵ���㣬��ֹѭ��
    if k==totalStep
       continue; 
    end
        
    % ���Ϳ������󣬼���ʵ�ʿ�����
    actualControlValue(k)=U_control(k-sigma_SC_reorder(k));
    
    X_state(:,k+1) = A_model_discrete * X_state(:,k) + B_model_discrete * U_control(k);
    Y_outputValue(k+1) = C_model_discrete * X_state(:,k+1);
    
end

% ʵ��ʹ�ÿ���������
figure(2)
stairs((1:totalStep)*T,U_control(1:totalStep),'-.k');
grid on
hold on

% ϵͳ�������
figure(3)
plot((1:totalStep)*T,Y_outputValue(1:totalStep),'-.k');
grid on
hold on
%}

% ********************************************************************************** %
% �� 1 ��
% ���⣺ʱ�Ӻ�����
% ��ʽ��ִ�������������򣨵�������ݰ���ѡȡʱ������µ����ݰ�������������������������������
for k=1:totalStep
    % ��ʼʱ�̣�����Ӧֵ���г�ʼ��
    if k==1
        % ��������ʼ��
        U_control=zeros(1,simulationTime);
        U_control(:,1)=0;
        % ϵͳ״ֵ̬��ʼ��
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % ϵͳ���ֵ��ʼ��
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(:,1)=C_model_discrete*X_state(:,1);
        % ����������ʼ��
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % ��ʼԤ��ֵ��ʼ��
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
        % ʵ��ʹ�ÿ���������ʼ��
        actualControlValue=zeros(1,simulationTime);
        % У��������ʼ��
        h=ones(N,1);
        % �������ʼ��
        error=zeros(1,simulationTime);
        % Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        Y_setValue=ones(P,1)*Sv;
        continue
    end
    
    % ���ݻ���������Ԥ�������У����������λ������ʱ�ӣ�������ʹ���յ����������ֵ����У��
    % ������ݰ����������������������˴�����
    if timeSequenceInt(k)==(d+1)
        % ������������ݰ�
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k);
        U_control(:,k)=U_control(:,k-1);
    else
        error(k)=Y_outputValue(k)-Y_predictedValue(1,k);
        tempCorrection=Y_predictedValue(:,k)+correction_h*h*error(k);
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % ������λ����S
        Y_predictedValue(:,k)=S*tempCorrection;
        % ���ݻ���������������������ݴ�����-��������ʱ�ӣ�ѡȡ���ʵ�A��������������
        controlIncrement(:,k)=K_gather(1+timeSequenceInt(k),:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        % ���ݻ������������������������ʱ��
        % controlIncrement(:,k)=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,k));
        
        % ���ݳ�ʼ����������������һʱ������ĳ�ʼԤ��ֵ
        Y_predictedValue(:,k+1)=Y_predictedValue(:,k)+stepResponse(:,1)*controlIncrement(1,k);
        
        U_control(:,k)=U_control(:,k-1)+controlIncrement(1,k);
    end
    
    % ĩβʱ�̣��������ʱ�̵����ֵ���㣬��ֹѭ��
    if k==totalStep
       continue; 
    end
        
    % ���Ϳ������󣬼���ʵ�ʿ�����
    actualControlValue(k)=U_control(k-sigma_SC_reorder(k));
    
    X_state(:,k+1) = A_model_discrete * X_state(:,k) + B_model_discrete * actualControlValue(k);
    Y_outputValue(k+1) = C_model_discrete * X_state(:,k+1);
    
end

% ʵ��ʹ�ÿ���������
figure(2)
stairs((1:totalStep)*T,U_control(1:totalStep),'-.k');
grid on
hold on

% ϵͳ�������
figure(3)
plot((1:totalStep)*T,Y_outputValue(1:totalStep),'-.k');
grid on
hold on

% ********************************************************************************** %
% �� 2 ��
% ���⣺ʱ�Ӻ�����
% ��ʽ��ִ���������������򣻿�������������������������
% ********************************************************************************** %
% �� 3 ��
% ���⣺ʱ�Ӻ�����
% ��ʽ��ִ�������������򣨵�������ݰ���ѡȡʱ������µ����ݰ�����������������������������������
% ********************************************************************************** %
% �� 4 ��
% ���⣺ʱ�Ӻ�����
% ��ʽ��ִ����������������ѡȡ���µ�������ݰ�����������������������������������
for k=1:totalStep
    % ��ʼʱ�̣�����Ӧֵ���г�ʼ��
    if k==1
        index=k+1;
        % ��������ʼ��
        U_control=zeros(1,simulationTime);
        U_control(:,1)=0;
        % ϵͳ״ֵ̬��ʼ��
        X_state=zeros(2,simulationTime);
        X_state(:,1)=[0 0]';
        % ϵͳ���ֵ��ʼ��
        Y_outputValue=zeros(1,simulationTime);
        Y_outputValue(:,1)=C_model_discrete*X_state(:,1);
        % ����������ʼ��
        controlIncrement=zeros(M,simulationTime);
        controlIncrement(:,1)=zeros(M,1);
        % ��ʼԤ��ֵ��ʼ��
        Y_predictedValue=zeros(N,simulationTime);
        Y_predictedValue(:,1)=Y_outputValue(1)*ones(N,1);
        % ʵ��ʹ�ÿ���������ʼ��
        actualControlValue=zeros(1,simulationTime);
        % У��������ʼ��
        h=ones(N,1);
        % �������ʼ��
        error=zeros(1,simulationTime);
        % Y_setValue=[Sv Sv ... Sv]',�ο��켣��������ֵ
        Y_setValue=ones(P,1)*Sv;
        continue
    end
    
    % �����ݰ�����Ϳ�ʼ����������
    while index<=totalStep && timeStampList_S2C(index)<=k
        error(index)=Y_outputValue(timeStampList_S2C(index))-Y_predictedValue(1,index);
        tempCorrection=Y_predictedValue(:,index)+correction_h*h*error(index);
        S=[zeros(N-1,1) eye(N-1);zeros(1,N-1) 1]; % ������λ����S
        Y_predictedValue(:,index)=S*tempCorrection;

        % ���ݻ������������������A���󲻸���ʱ�Ӷ�̬�仯
        controlIncrement(:,timeStampList_S2C(index))=K_gather(1,:)*(Y_setValue-Y_predictedValue(1:P,index));

        U_control(:,index)=U_control(:,index-1)+controlIncrement(1,index);
        % ���ݳ�ʼ����������������һʱ������ĳ�ʼԤ��ֵ
        Y_predictedValue(:,index+1)=Y_predictedValue(:,index)+stepResponse(:,1)*controlIncrement(1,index);
        
        index=index+1;
    end
    
    % ĩβʱ�̣��������ʱ�̵����ֵ���㣬��ֹѭ��
    if k==totalStep
       continue; 
    end
    
    % ���Ϳ������󣬼��������ʵ��ʹ�õĿ�����
    actualControlValue(k)=U_control(k-latestArrive_SC_noreorder(k));
    
    X_state(:,k+1) = A_model_discrete * X_state(:,k) + B_model_discrete * actualControlValue(k);
    Y_outputValue(k+1) = C_model_discrete * X_state(:,k+1);
    
end

% ʵ��ʹ�ÿ���������
figure(2)
set(gcf,'Position',[1000,400,700,500]);
stairs((1:totalStep)*T,U_control(1:totalStep),'--k');
legend('��׼���','��������','����������');
xlabel('ʱ��k');
ylabel('ʵ��ʹ�ÿ�����u');
grid on
hold on

% ϵͳ�������
figure(3)
set(gcf,'Position',[1000,400,700,500]);
plot((1:totalStep)*T,Y_outputValue(1:totalStep),'--k');
legend('��׼���','��������','����������');
xlabel('ʱ��k');
ylabel('���y');
grid on
hold on

