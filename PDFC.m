function [ output ] = k_means(data,cen,quan,means, k_value)
% ���ܣ�ʵ��PDFC�㷨�ľ��๦�ܣ�
% ���룺    data, Ϊһ�� ���� M��N�� ��ʾ������������M��ʾ����M����������N��ʾÿһ��������ά�ȣ�
%           k_value, ��ʾ����������Ŀ��
% �����    output, ��һ�������� M��������ʾÿһ���������ڵ�����ţ�


%�������У����ѡȡK��������Ϊ��ʼ�ľ������ģ�
data_num = size(data, 1);
center = cen;
%���ڼ�������������
iteration = 0;
while 1
    %�����������������ĵľ��룻
    distance = euclidean_distance(data, center);
    %����������ÿһ�д�С�������� �����Ӧ��indexֵ����ʵ����ֻ��Ҫindex�ĵ�һ�е�ֵ��
    [~, index] = sort(distance, 2, 'ascend');
    for i=1:length(quan)
        for j=1:k_value
          distance = euclidean_distance(data(i,:), center(j,:));
          F(i,j)=(quan(i)*means(j))/distance.^2;%��ÿ�����ĵ㵽��ֵ��������С
        end
    end
    [~, index] = sort(F, 2, 'descend');
    %�������γ��µľ������ģ�
    xc = index(:, 1);
     xvc=cell(1,k_value)
    for i=1:k_value
        for j=1:length(quan)
        if xc(j)==i  
        xvc{i}=[xvc{i};quan(j)];
        end
        end
    end
    
    for i=1:k_value
         meanx(i)=mean(xvc{i});
    end
    center_new = zeros(k_value, size(data, 2));
    for i = 1:k_value
        data_for_one_class = data(index(:, 1) == i, :);          
        center_new(i,:) = mean(data_for_one_class, 1);    %��Ϊ��ʼ�ľ�������Ϊ�������е�Ԫ�أ����Բ������ĳ������������Ϊ0�������
    end
   
   
    
    %����������������۾�һ��������
    iteration = iteration + 1;
    fprintf('The number of iterations is��%d\n', iteration);
    dsq=sum(abs(means-meanx));
    % ��������εľ������Ĳ��䣬��ֹͣ����������ѭ����
    sd=sum(sum(abs(center_new-center)));
    fprintf('%d\n', sd);
      fprintf('%d\n', dsq);
    if (center_new == center)&(dsq<0.01)
        break;
    end
    means=meanx;
    center = center_new;
    
    
 
    
    
end

output = index(:, 1);
    
end