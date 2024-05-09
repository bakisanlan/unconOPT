function x_new = line_search(x, d)

% finding next x variable with using linear search method
% x_new, x and d is Nx1 column vector

syms a
x_new = x + a .* d;
end