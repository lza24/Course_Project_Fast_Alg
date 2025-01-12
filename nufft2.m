function f = nufft2(c, x, y, w, n)
% 二维 nufft 求解器, 采用 Kaiser-Bessel 作为 spreading function
% x,y 为采样点的坐标, 均在 [-pi, pi] 内
% c 为 x,y 处对应的采样值
% w 为 spreading function 的支集宽度, 表示支集内所含过采样点的个数
% n 为每个方向上 DFT 规模
% 输出为 [-n/2,n/2]*[-n/2,n/2] 上的 DFT

M = length(x);        % 读取采样点个数
h = 2 * pi / n;       % 过采样间距
alpha = pi * w / n;   % 支集范围 [-alpha, alpha]
beta = 2.3 * w;       % spreading function 的参数, 表示截断的程度, 一般设置为 2.3w
b = zeros(n,n);       % 初始化过采样网格
f = zeros(n,n);       % 初始化输出结果

% 利用 spreading function 作卷积将非均匀采样点过渡到均匀的过采样网格上
for l_1 = 0:1:n-1
    for l_2 = 0:1:n-1
        for j = 1:1:M 
            b(l_1+1, l_2+1) = b(l_1+1, l_2+1) + c(j) * phi2(alpha, beta, l_1*h - x(j), l_2*h - y(j));
        end
    end
end

% 在过采样网格上使用 fft
b_hat = fft2(b);

% 系数修正
for k_1 = 0:1:n-1
    for k_2 = 0:1:n-1
        phi_hat1 = 2 * sinh(sqrt(beta^2 - (-n/2 + k_1)^2 * alpha^2)) / (besseli(0, beta) * sqrt(beta^2 - (-n/2 + k_1)^2 * alpha^2));
        phi_hat2 = 2 * sinh(sqrt(beta^2 - (-n/2 + k_2)^2 * alpha^2)) / (besseli(0, beta) * sqrt(beta^2 - (-n/2 + k_2)^2 * alpha^2));
        p = 4 / (w^2 * phi_hat1 * phi_hat2);
        f(k_1+1, k_2+1) = p * b_hat(k_1+1, k_2+1);
    end
end

% 移位频谱，使零频率成分位于中心, 符合正常习惯
f = fftshift(f);

end