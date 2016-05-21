% 输出结果为矩阵stepResponse，维数为N*(d+1)，每列依次存储着时延为0:d是系统对应的阶跃系数构成的列向量

clear all;clc;close all;

T=0.04; % 采样周期
N=40; % 模型长度（建模时域）
d=5; % 最大时延周期数

% 系统的连续状态空间方程
A_model=[-4 -0.03;0.75 -10];
B_model=[2 0]';
C_model=[0 1];

% 系统的离散状态空间方程
[A_model_discrete,B_model_discrete]=c2d(A_model,B_model,T);
C_model_discrete=C_model;

% 无时延系统的阶跃响应，获得阶跃系数
temp_response=dstep(A_model_discrete,B_model_discrete,C_model_discrete,0);
stepResponse(:,1)=[zeros(N,1),eye(N),zeros(N,length(temp_response)-(N+1))]*temp_response;

% 时延为delay=1:d的系统的阶跃系数
for delay=1:d 
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