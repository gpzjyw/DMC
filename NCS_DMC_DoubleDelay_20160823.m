%  2016-08-23�����������¸Ľ���1.����˫��ʱ�ӵ�������Ҷ����ܳ�����������⣻
%                             2.����Ĵ���ִ����ѡȡ���¿����źţ�������Ҳѡ�����µ����ֵ��
%                             3.��������У�����֣����õĿ�����������ʱ�� ������-������ ʱ��ȷ����
%                             4.�����������¼��������յ�һ�����ֵ����һ���¼�����ѡ�õ����ֵΪ���µ����ֵ��
% ********************************************************************************** %

clear all;clc;close all;

d=4; % ���ʱ��Ϊ d*T
T=0.04; % ��������
totalStep=50;% �����ܲ���
timeSequence = d*T*[0 rand(1,totalStep-d-1)-0.2 zeros(1,d)]; % ÿ��ʱ�����ݰ���ʱ�Ӵ�С������ʱ�ӵĸ���Ϊ80%��
timeSequenceInt = timeSequence; % �������Ĵ��ڣ���ʱ��ת��Ϊ�̶���������ʱ�ӣ�0,T,2T,...,dT
timeSequenceLength = length(timeSequence);
controlSequence=zeros(2*(d+1),totalStep);

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

% �趨�����ĸ��ʣ�������������ʱ�̴���dropoutPoint����������ʱ�̵�ʱ������Ϊ(d+1)T�����������ĸ���Ϊ10%
dropoutRate=[0 rand(1,totalStep-d-1) zeros(1,d)];
dropoutPoint=[]; % ����ʱ�̵ļ���
for i=1:totalStep
    if dropoutRate(i)>=0.9
        timeSequenceInt(i)=d+1;
        dropoutPoint=[dropoutPoint,i];
    end
end

% ȷ��ÿ��ʱ��֮ǰ������������
continuousDropoutNum=zeros(1,totalStep);
for i=1:length(dropoutPoint)
    if continuousDropoutNum(dropoutPoint(i))~=0
        continuousDropoutNum(dropoutPoint(i)+1)=continuousDropoutNum(dropoutPoint(i))+1;
    else
        continuousDropoutNum(dropoutPoint(i)+1)=1;
    end
end

% ����ʱ��ȷ������ʱ�̵ļ���
disorderingPoint=[]; % ����ʱ�̵ļ���
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

%ȡ��δ����������ʱ�̵�ʱ�ӣ�����timeSequenceIntDraw�У�����δ������ʱ�̴���delayPoint��
timeSequenceIntDraw=[];
j=1;
delayPoint=[]; % ʱ��ʱ�̵ļ���
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

% ����sigma(k)��ִ�������ݸ�ֵѡȡ���µĿ�����
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

% ��ͼ
figure(1)
t=T*(1:totalStep);
subplot(2,1,1);
plot(delayPoint*T,timeSequenceIntDraw*T,'*k',disorderingPoint*T,zeros(length(disorderingPoint),1),'ok',dropoutPoint*T,d*T,'xk');
axis([0 totalStep*T -0.5*T (d+1)*T]);
legend('ʱ��','���������ʱ��','����������');
xlabel('ʱ��k');
ylabel('\tau _k');
grid on
hold on
subplot(2,1,2);
plot((1:totalStep)*T,sigma,'*k')
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




























