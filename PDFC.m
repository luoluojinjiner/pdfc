function [ output ] = k_means(data,cen,quan,means, k_value)
% 功能：实现PDFC算法的聚类功能；
% 输入：    data, 为一个 矩阵 M×N， 表示样本集，其中M表示共有M个样本，　N表示每一个样本的维度；
%           k_value, 表示聚类的类别数目；
% 输出：    output, 是一个列向量 M×１，表示每一个样本属于的类别编号；


%从样本中，随机选取K个样本作为初始的聚类中心；
data_num = size(data, 1);
center = cen;
%用于计数迭代次数：
iteration = 0;
while 1
    %获得样本集与聚类中心的距离；
    distance = euclidean_distance(data, center);
    %将距离矩阵的每一行从小到大排序， 获得相应的index值，其实我们只需要index的第一列的值；
    [~, index] = sort(distance, 2, 'ascend');
    for i=1:length(quan)
        for j=1:k_value
          distance = euclidean_distance(data(i,:), center(j,:));
          F(i,j)=(quan(i)*means(j))/distance.^2;%算每个核心点到峰值的引力大小
        end
    end
    [~, index] = sort(F, 2, 'descend');
    %接下来形成新的聚类中心；
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
        center_new(i,:) = mean(data_for_one_class, 1);    %因为初始的聚类中心为样本集中的元素，所以不会出现某类别的样本个数为0的情况；
    end
   
   
    
    %输出迭代次数，给眼睛一个反馈；
    iteration = iteration + 1;
    fprintf('The number of iterations is：%d\n', iteration);
    dsq=sum(abs(means-meanx));
    % 如果这两次的聚类中心不变，则停止迭代，跳出循环；
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