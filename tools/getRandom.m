function result = getRandom( startNum, endNum )
% 随机生成给定范围内的数字
% 给定起始和结束，随机生成该范围内的数字（不一定为整数）

result = startNum + rand(1)*(endNum-startNum);

end