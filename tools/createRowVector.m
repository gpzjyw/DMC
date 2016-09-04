function [ arr ] = createRowVector ( index, M )
% 根据输入值 index ，生成行向量，要保证第 (index-1) 为 1 ，其余项为 0

arr=zeros(1,M);
arr(index+1)=1;

end

