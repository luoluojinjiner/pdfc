function [ output ] = euclidean_distance(data, center)
% ���ڼ���ѵ��������������ĵĵ�ŷ�Ͼ����ƽ����
% ����  dataΪһ�� ���� M��N�� ��ʾ������������M��ʾ����M����������N��ʾÿһ��������ά�ȣ�
%      centre Ϊһ������ K��N����ʾK���������ģ�N��ʾ������ά�ȣ�
%      output Ϊһ�����󣬴�СΪM��K�� �ڣ��У��б�ʾ��X���������Y���������ĵľ��룻��ÿһ�б�ʾһ��������K���������ĵľ��룩��


% ���ߣ�����壻
% ʱ�䣺2017��10��14�գ�




data_num = size(data, 1);
center_num = size(center, 1);
output = zeros(data_num, center_num);
for i = 1:center_num
    difference = data - repmat(center(i,:), data_num, 1);    %�����������i���������ĵĲ
    sum_of_squares = sum(difference .* difference, 2);        %��ƽ���� ����ÿһ����ͣ�
    output(:, i) = sum_of_squares;             
end

end