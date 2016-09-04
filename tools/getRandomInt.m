function result = getRandomInt( startNum, endNum )
% 随机生成给定范围内的整数
% 给定起始整数和结束整数，随机生成该范围内的整数

result = startNum + round(rand(1)*(endNum-startNum));

end