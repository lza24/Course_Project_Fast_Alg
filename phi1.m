function phi = phi1(alpha, beta, x, y)
% Guassian 函数, 用于 spreading function
% x, y 为输入的两个分量
% 截断范围为 [-alpha, alpha]
% 截断处的误差为 exp(-beta)

% 周期延拓
x = mod(x, 2 * pi);
y = mod(y, 2 * pi);
if abs(x) <= alpha && abs(y) <= alpha
    phi = exp(-beta * (x.^2 + y.^2) / alpha^2);
else
    phi = 0;

end

