function [ arr ] = createRowVector ( index, M )
% ��������ֵ index ��������������Ҫ��֤�� (index-1) Ϊ 1 ��������Ϊ 0

arr=zeros(1,M);
arr(index+1)=1;

end

