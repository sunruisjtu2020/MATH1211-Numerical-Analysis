% This Program is distributed under the terms of MIT license.
% Author: Rui Sun
% Date: 2021-12-12
% Copyright: Rui Sun from SAA, SJTU
% Description: Task 1 - Lagrange-Polynomial
% Environment: MATLAB R2021b on Windows

% main
n = [5, 10, 20];
x = [-0.95, -0.47, 0.1, 0.37, 0.93];

% function f(x) = 1 / (1 + 25x^2)
for i = n
    fprintf('n = %d\n', n);
    for j = x
        fprintf('f(%f)-p_n(%f)=%f\n', j, j, f(j) - L(i, j, f));
    end
end

for i = n
    fprintf('n = %d\n', n);
    for j = x
        fprintf('g(%f)-p_n(%f)=%f\n', j, j, g(j) - L(i, j, g));
    end
end

for i = n
    fprintf('n = %d\n', n);
    for j = x
        fprintf('f(%f)-p_n(%f)=%f\n', j, j, f(j) - redefinedL(i, j, f));
    end
end


% Lagrange-Polynomial
function res = L(n, x, f)
    x_ = linspace(-1, 1, n + 1);
    y_ = f(x_);
    l_ = zeros(size(x_, 2));
    res = 0;
    for i = 1:n + 1
        deno = 1;
        nume = 1;
        for j = 1:n + 1
            if i ~= j
                deno = deno * (x_(i) - x_(j));
                nume = nume * (x - x_(j));
            end
        end
        l_(i) = nume / deno;
        res = res + l_(i) * y_(i);
    end
end

function res = LLL(n, x, f)
    res = zeros(200);
    for i = 1:200
        res(i) = L(n, x(i), f);
    end
end

function res = redefinedL(n, x, f)
    x_ = zeros(n);
    for i = 1:n
        x_(i) = cos((2 * i + 1) * pi / (2 * (n + 1)));
    end
    y_ = f(x_);
    l_ = zeros(size(x_, 2));
    res = 0;
    for i = 1:n + 1
        deno = 1;
        nume = 1;
        for j = 1:n + 1
            if i ~= j
                deno = deno * (x_(i) - x_(j));
                nume = nume * (x - x_(j));
            end
        end
        l_(i) = nume / deno;
        res = res + l_(i) * y_(i);
    end
end

function res = redefinedLLL(n, x, f)
    res = zeros(200);
    for i = 1:200
        res(i) = redefinedL(n, x, f);
    end
end

function res = f(x)
    res = 1 / (1 + 25 * x^2);
end

function res = g(x)
    res = exp(x);
end