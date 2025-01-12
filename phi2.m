function phi = phi2(alpha, beta, x, y)
% Kaiser-Bessel 函数, 用于 spreading function
% x, y 为输入的两个分量
% 截断范围为 [-alpha, alpha]
% beta 为函数参数

% 周期延拓
x = mod(x, 2 * pi);
y = mod(y, 2 * pi);
if abs(x) <= alpha && abs(y) <= alpha
    phi = besseli(0, beta * sqrt(1 - x.^2 / alpha^2)) ...
        * besseli(0, beta * sqrt(1 - y.^2 / alpha^2)) / besseli(0, beta)^2;
else
    phi = 0;

end
