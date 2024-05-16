clc
close all
clear all

data = xlsread("data.xlsx");

%专家权重
e_j = xlsread("data1.xlsx");

seita = 36 %1,2.25,10,36

Z_1 = data(1:5,:);
Z_2 = data(7:11,:);
Z_3 = data(13:end,:);

% Wjk 公式6 最大化 
sum_1 = sum(Z_1(:));
sum_2 = sum(Z_2(:));
sum_3 = sum(Z_3(:));
for i =1:5
    W_j_k1(i) = sum(Z_1(:,i)) / sum_1;
end
for i =1:5
    W_j_k2(i) = sum(Z_2(:,i)) / sum_2;
end
for i =1:5
    W_j_k3(i) = sum(Z_3(:,i)) / sum_3;
end

% Wj*k 公式7的值，用的就是这个wj*k，三个专家分别对应有omega_j_1，omega_j_2，omega_j_3
n = 5; 
m = 5;
% 初始化分子和分母
numerator = 0;
denominator = 0;
% 计算分子和分母
for j = 1:n
    for i = 1:m
        for h = 1:m
            numerator = numerator + abs(Z_1(i, j) - Z_1(h, j));
        end
    end
    numerator_1(j) = numerator;
    numerator = 0;
    for i = 1:m
        for h = 1:m
            denominator = denominator + abs(Z_1(i, j) - Z_1(h, j));
        end
    end
end
% 计算 omega_j_1
omega_j_1 = numerator_1 ./ denominator;

% 初始化分子和分母
numerator = 0;
denominator = 0;
% 计算分子和分母
for j = 1:n
    for i = 1:m
        for h = 1:m
            numerator = numerator + abs(Z_2(i, j) - Z_2(h, j));
        end
    end
    numerator_1(j) = numerator;
    numerator = 0;
    for i = 1:m
        for h = 1:m
            denominator = denominator + abs(Z_2(i, j) - Z_2(h, j));
        end
    end
end
% 计算 omega_j_2
omega_j_2 = numerator_1 ./ denominator;

% 初始化分子和分母
numerator = 0;
denominator = 0;
% 计算分子和分母
for j = 1:n
    for i = 1:m
        for h = 1:m
            numerator = numerator + abs(Z_3(i, j) - Z_3(h, j));
        end
    end
    numerator_1(j) = numerator;
    numerator = 0;
    for i = 1:m
        for h = 1:m
            denominator = denominator + abs(Z_3(i, j) - Z_3(h, j));
        end
    end
end
% 计算 omega_j_3
omega_j_3 = numerator_1 ./ denominator;

% 计算任意两方案优势度，得三个专家的对应值
phi_k_j_1 = zeros(n, n); % 创建一个零矩阵来存储结果
for k = 1:n
    for t = 1:m
        for i = 1:m
            d = (Z_1(t,k) - Z_1(i,k))^2;
            if Z_1(i, k) > Z_1(t, k)
                phi_k_j_1(i, t) = phi_k_j_1(i, t) + sqrt((omega_j_1(1,k) * d) / sum(omega_j_1(:)));
            elseif Z_1(i, k) == Z_1(t, k)
                phi_k_j_1(i, t) = phi_k_j_1(i, t) + 0;
            else
                phi_k_j_1(i, t) = phi_k_j_1(i, t) - (1 / seita) * sqrt(( sum(omega_j_1(:))* d) / omega_j_1(1,k));
            end
        end
    end
    phi_1{k} = phi_k_j_1;
end

phi_k_j_2 = zeros(n, n); % 创建一个零矩阵来存储结果
for k = 1:n
    for t = 1:m
        for i = 1:m
            d = (Z_1(t,k) - Z_1(i,k))^2;
            if Z_1(i, k) > Z_1(t, k)
                phi_k_j_2(i, t) = phi_k_j_2(i, t) + sqrt((omega_j_2(1,k) * d) / sum(omega_j_2(:)));
            elseif Z_1(i, k) == Z_1(t, k)
                phi_k_j_2(i, t) = phi_k_j_2(i, t) + 0;
            else
                phi_k_j_2(i, t) = phi_k_j_2(i, t) - (1 / seita) * sqrt(( sum(omega_j_1(:))* d) / omega_j_1(1,k));
            end
        end
    end
    phi_2{k} = phi_k_j_2;
end

phi_k_j_3 = zeros(n, n); % 创建一个零矩阵来存储结果
for k = 1:n
    for t = 1:m
        for i = 1:m
            d = (Z_1(t,k) - Z_1(i,k))^2;
            if Z_1(i, k) > Z_1(t, k)
                phi_k_j_3(i, t) = phi_k_j_3(i, t) + sqrt((omega_j_3(1,k) * d) / sum(omega_j_3(:)));
            elseif Z_1(i, k) == Z_1(t, k)
                phi_k_j_3(i, t) = phi_k_j_3(i, t) + 0;
            else
                phi_k_j_3(i, t) = phi_k_j_3(i, t) - (1 / seita) * sqrt(( sum(omega_j_1(:))* d) / omega_j_1(1,k));
            end
        end
    end
    phi_3{k} = phi_k_j_3;
end

% 计算总体优势度 
result_matrix = zeros(5, 5);
% 遍历 phi_1 中的每个 cell
for k = 1:5 
    % 将当前 cell 中的数组与结果矩阵相加
    result_matrix = result_matrix + phi_1{1, k};
end

% 计算总体优势度 第一个公式的第一个专家的值
phi_j_1 = result_matrix

result_matrix = zeros(5, 5);
% 遍历 phi_2 中的每个 cell
for k = 1:5 
    result_matrix = result_matrix + phi_2{1, k};
end
% 计算总体优势度 第一个公式的第二个专家的值
phi_j_2 = result_matrix

result_matrix = zeros(5, 5);
% 遍历 phi_3 中的每个 cell
for k = 1:5 
    result_matrix = result_matrix + phi_3{1, k};
end
% 计算总体优势度 第一个公式的第三个专家的值
phi_j_3 = result_matrix

% 计算总体优势度第二个公式
phi_it_1 = phi_j_1 .* e_j(1,1);
phi_it_2 = phi_j_2 .* e_j(1,2);
phi_it_3 = phi_j_3 .* e_j(1,3);
% 计算总体优势度 第二个公式的的值
phi_it = phi_it_1 + phi_it_2 + phi_it_3

% 计算总体优势度 第三个公式
phi_Z_i = sum(phi_it,2);
phi_Z = (phi_Z_i - min(phi_Z_i)) / (max(phi_Z_i) - min(phi_Z_i))



