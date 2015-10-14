function  drawDMCPicture = DMC_func(Ts,N,P,M,Sv,control_R,error_Q,correction_h,y0)
    % Ts，采样时间
    % N，模型长度
    % P，预测时域
    % M，控制时域
    % Sv，设定值，参考轨迹
    % control_R，控制权矩阵系数
    % error_Q，误差权矩阵系数
    % correction_h，校正系数 
    % y0,阶跃信号输入

    %控制矩阵A、A0,控制矩阵K由下段程序得出：
    A=zeros(P,M);
    a=zeros(N,1);
    for i=1:N
        a(i)=y0(i+1);%y0为系统阶跃响应的系数
    end
    for i=1:P
        for j=1:M
            if i-j+1>0
                A(i,j)=a(i-j+1);
            end
        end
    end

    K=inv(A'*eye(P)*error_Q*A+eye(M)*control_R)*A';
    ys=ones(N,1);
    y=zeros(N,1);
    u=zeros(N,1);
    e=zeros(N,1);
    A0=zeros(P,N-1);
    %A0为系统模型，用于计算准确的输出值
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

    %循环部分，变量初始值为0
    %DMC程序
    for k=2:N
        Uk_1=zeros(N-1,1);
        for i=1:(N-1)
            if k-N+i<=0
                Uk_1(i)=0;
            else
                Uk_1(i)=u(k-N+i);
            end
        end
        Y0=A0*Uk_1;
        e(k)=y(k-1)-Y0(1);%计算误差
        Ysk=zeros(P,1);
        for i=1:P
            Ysk(i)=Sv;%Ysk=[Sv Sv ... Sv]'
        end
        Ek=zeros(P,1);
        for i=1:P
            Ek(i)=e(k);%Ek=[e(k) e(k) ...e(k)]'
        end
        dersu=K*(Ysk-Y0-Ek);%计算控制增量
        for i=1:M
            if k+i-1<=N
                u(k+i-1)=u(k+i-1-1)+dersu(i);%控制量更新，滚动优化
            end
        end
        temp=0;
        for j=1:N-1
            if k-j<=0
                temp;
            else
                if k-j-1<=0  
                    temp=temp+a(j)*u(k-j);
                else
                    temp=temp+a(j)*(u(k-j)-u(k-j-1));%预测输出
                end
            end
        end
        if k-N<=0
            y(k)=temp+e(N)*correction_h;%基于误差的预测输出校正
        else
            y(k)=temp+a(N)*u(k-N)+e(N)*correction_h;%基于误差的预测输出校正
        end
    end

    %画图显示结果
    figure(2)
    t=Ts.*(1:N);
    subplot(2,1,1);
    plot(t,y,'.-');
    legend('y','Location','Best');
    title('输出曲线');
    xlabel('Time')
    ylabel('振幅 ')
    grid on
    subplot(2,1,2);
    plot(t,u,'.-');
    legend('控制作用u','Location','Best');
    title('控制作用');
    xlabel('Time')
    ylabel('振幅')
    grid on
    hold on

end