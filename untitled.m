clc
close all
clear all

%相对权重
data = xlsread("data.xlsx");

Z_1 = data(1:5,:);
Z_2 = data(7:11,:);
Z_3 = data(13:end,:);

%第一步
sum_1 = sum(Z_1(:));
sum_2 = sum(Z_2(:));
sum_3 = sum(Z_3(:));

b_1 = Z_1;
b_2 = Z_3;
b_3 = Z_2;

%第二步
numerator_1 = nchoosek(3, 1);
numerator_2 = nchoosek(3, 2);
numerator_3 = nchoosek(3, 3);
sum_1 = numerator_1 + numerator_2 + numerator_3;
beita_1 = numerator_1 / sum_1;
beita_2 = numerator_2 / sum_1;
beita_3 = numerator_3 / sum_1;

%第三步
w_1 = beita_1 .* b_1;
w_2 = beita_2 .* b_2;
w_3 = beita_3 .* b_3;
w_r = [w_1,w_2,w_3];
w_r_1 = sum(w_r,2);

%第四步
for i=1:5
    ur(i,1) = w_r_1(i,1) / sum(w_r_1);
end
















