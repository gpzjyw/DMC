function [ x_state_next,y_output_next ] = stateSpace( x_state_k,u_control_k,A,B,C )
% ����״̬�ռ䷽�̼������
% ����kʱ�̵�״ֵ̬�Ϳ�����������k+1ʱ�̵�״ֵ̬�����ֵ

x_state_next=A*x_state_k+B*u_control_k;
y_output_next=C*x_state_k;

end

