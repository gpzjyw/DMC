function result = getRandomInt( startNum, endNum )
% ������ɸ�����Χ�ڵ�����
% ������ʼ�����ͽ���������������ɸ÷�Χ�ڵ�����

result = startNum + round(rand(1)*(endNum-startNum));

end