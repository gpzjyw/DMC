% ������Ϊ����stepResponse��ά��ΪN*(d+1)��ÿ�����δ洢��ʱ��Ϊ0:d��ϵͳ��Ӧ�Ľ�Ծϵ�����ɵ�������

clear all;clc;close all;

T=0.04; % ��������
N=40; % ģ�ͳ��ȣ���ģʱ��
d=5; % ���ʱ��������

% ϵͳ������״̬�ռ䷽��
A_model=[-4 -0.03;0.75 -10];
B_model=[2 0]';
C_model=[0 1];

% ϵͳ����ɢ״̬�ռ䷽��
[A_model_discrete,B_model_discrete]=c2d(A_model,B_model,T);
C_model_discrete=C_model;

% ��ʱ��ϵͳ�Ľ�Ծ��Ӧ����ý�Ծϵ��
temp_response=dstep(A_model_discrete,B_model_discrete,C_model_discrete,0);
stepResponse(:,1)=[zeros(N,1),eye(N),zeros(N,length(temp_response)-(N+1))]*temp_response;

% ʱ��Ϊdelay=1:d��ϵͳ�Ľ�Ծϵ��
for delay=1:d 
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