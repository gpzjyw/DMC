A=[0.91 1 0.5 0.5 ; 0 0.91 1 1 ; 0 0 0.91 0 ; 0 0 0 0.606];
B=[0 0 0 0.00792]';
C=[1 0.5 0 0];
%计算时间t=0-200内的状态值，u=0
for(t=0:1:200)
    if(t==0)
        temp_x1=[1 1 1 1]';
        temp_x=temp_x1;
    else
        temp_x2=A*temp_x1;
        temp_x1=temp_x2;
        temp_x=[temp_x temp_x2];
    end    
end
%计算时间t=0-200内的输出值
t=0:1:200;
temp_y=C*temp_x;

subplot(2,1,1);
plot(t,temp_y);
title('Y');

subplot(2,1,2);
plot(t,temp_x(1,:),t,temp_x(2,:),t,temp_x(3,:),t,temp_x(4,:));
title('X');