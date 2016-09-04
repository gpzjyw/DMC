function [ x_state_next,y_output_next ] = stateSpace( x_state_k,u_control_k,A,B,C )
% 根据状态空间方程计算输出
% 基于k时刻的状态值和控制量，计算k+1时刻的状态值和输出值

x_state_next=A*x_state_k+B*u_control_k;
y_output_next=C*x_state_k;

end

