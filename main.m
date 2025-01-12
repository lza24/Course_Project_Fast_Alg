%% 偏微分方程快速算法 大作业
% 非均匀快速傅里叶变换 (NUFFT)
% 刘正傲 2024312443

%
%% spreading function 的绘制
% 选用两种不同的 spreading function, Guassian 和 Kaiser-Bessel
% 绘制他们的示意图
alpha = 1;
beta = 5;
x = -1:0.01:1;
phi_1 = exp(-0.5 * beta * (x.^2 / alpha^2));
phi_2 = besseli(0, beta * sqrt(1 - x.^2 / alpha^2)) / besseli(0, beta);
figure;
plot(x, phi_1, 'LineWidth', 2);
xlabel('x');
ylabel('\phi(x)');
title('Spreading function: Guassian');
grid on;
figure;
plot(x, phi_2, 'LineWidth', 2);
xlabel('x');
ylabel('\phi(x)');
title('Spreading function: Kaiser-Bessel');
grid on;
%}


%
%% 使用不同 spreading function 的性能测试(计算复杂度验证)
% 针对不同的 DFT 数量 N, 进行多次试验, 分别记录运行时间

time1 = zeros(1,10);
time2 = zeros(1,10);
N = (4 * 2.^(1:1:10)).^2;      % DTF 的规模
for i = 1:1:10
    % 随机生成采样点
    x = -pi + 2 * pi * rand(10, 1);  % x 的范围是 [-pi, pi]
    y = -pi + 2 * pi * rand(10, 1);  % y 的范围是 [-pi, pi]

    % 生成采样值
    c = sin(3 * x) .* sin(1 * y);

    % 使用 Guassian 作为扩展函数的求解器计算并计时
    tic;
    w = 10;     
    f1 = nufft1(c, x, y, w, sqrt(N(i)));
    time1(i) = toc;

    % 使用 Kaiser-Bessel 作为扩展函数的求解器计算并计时
    tic;
    w = 10;     
    f2 = nufft2(c, x, y, w, sqrt(N(i)));
    time2(i) = toc;
end

% 作图
figure;
loglog(N, time1, 'LineWidth', 2);
xlabel('Number of DFT points');
ylabel('Execution time');
title('Execution time with number of DFT points');

figure;
loglog(N, time2, 'LineWidth', 2);
xlabel('Number of DFT points');
ylabel('Execution time');
title('Execution time with number of  DFT points');
%}


%
%% 实例应用(信号处理)
% 以 f = sin(2x)sin(y) + sin(4x)sin(3y) 为例, 展示结果

M = 2000;  % 随机样点个数
n = 64;    % 每个维度 DFT 的个数

% 随机生成 M 个点 (x, y) 作为采样点
x = -pi + 2 * pi * rand(M, 1);  
y = -pi + 2 * pi * rand(M, 1);  

% 生成采样值
c = sin(2 * x) .* sin(1 * y) + sin(4 * x) .* sin(3 * y);

% 计算其 DFT
w = 10; 
f = nufft1(c,x,y,w,n);

%作图
k_x = -n / 2 + (0:1:n-1); 
k_y = -n / 2 + (0:1:n-1);  
figure;
imagesc(k_x, k_y, abs(f));
colorbar;
title('frequency spectrum');
xlabel('frequency x');
ylabel('frequency y');
axis equal;

% 展示原始函数图像
x = linspace(-pi, pi, 50);  
y = linspace(-pi, pi, 50);  
[X, Y] = meshgrid(x, y);
f = sin(2 * X) .* sin(Y) + sin(4 * X) .* sin(3 * Y);
figure;
surf(X, Y, f);
title('f(x,y)');
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
colorbar;  

%}


