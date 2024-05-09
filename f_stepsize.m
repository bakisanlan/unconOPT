function f_stePsize = f_stepsize(f, x_new)
% this function convert any f(x1,x2,x3..) function to f(step_size).

% f is function handle, x_new is sym var

% x_new = num2cell(x_new);
% f_stePsize = @(a) f(x_new{:});
% 
% f_stePsize = f_stePsize(1);

f_stePsize = @(a) f(x_new);
f_stePsize = f_stePsize(1);

end