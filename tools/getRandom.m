function result = getRandom( startNum, endNum )
% ������ɸ�����Χ�ڵ�����
% ������ʼ�ͽ�����������ɸ÷�Χ�ڵ����֣���һ��Ϊ������

result = startNum + rand(1)*(endNum-startNum);

end