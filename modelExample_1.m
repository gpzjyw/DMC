%连续的状态空间方程
A=[-4 -0.03;0.75 -10];
B=[2 0]';
C=[0 1];

T=0.04;%采样周期
N=50;%模型时域

t=[0:T:N*T];
y0=step(ss(A,B,C,0),t);
figure(1)
plot(t,y0)
hold on

DMC_func(T ,N ,15 ,5 ,1 ,0.01 ,1 ,1,y0);
    % Ts，采样时间
    % N，模型长度
    % P，预测时域
    % M，控制时域
    % Sv，设定值，参考轨迹
    % control_R，控制权矩阵系数
    % error_Q，误差权矩阵系数
    % correction_h，校正系数 
    % y0,阶跃信号输入