% 指定你的CSV文件的路径
file_path = '/Users/zhouyaya/Desktop/犹豫度得分.csv';

% 使用csvread函数读取CSV文件
% 请根据CSV文件的实际情况来调整参数，例如，如果有标题行，可以使用csvread(file_path, 1, 0)来跳过标题行。
data = csvread(file_path);

% 现在，"data" 变量将包含你的Excel数据作为一个矩阵
totalSum = sum(data); % 对整个矩阵进行求和

disp(totalSum);

