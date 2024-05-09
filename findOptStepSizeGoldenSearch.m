function optStepSize = findOptStepSizeGoldenSearch(f, x_init, d)
% this function find optimal step size to minimize f(x_init) in the
% direction d.

%eps = 0.0001;
eps = 0.0005;
sigma = 0.005;
GR = (sqrt(5) - 1)/2;

x_new = line_search(x_init,d);
a = symvar(x_new);

f_stepS = f_stepsize(f,x_new);

%% phase 1 loop
f_alpha_k = 0;
f_alpha_k1 = 0;
f_alpha_k2 = 0;
q=0;
i = 0;
termination = 10;
while ~((f_alpha_k1 < f_alpha_k) && (f_alpha_k1 < f_alpha_k2))
    alpha_k  =  sumalpha(q,sigma);
    alpha_k1 = sumalpha(q+1,sigma);
    alpha_k2 = sumalpha(q+2,sigma);

    f_alpha_k  =  double(subs(f_stepS,a,alpha_k));
    f_alpha_k1 = double(subs(f_stepS,a,alpha_k1));
    f_alpha_k2 = double(subs(f_stepS,a,alpha_k2));

    q = q +1;

    % sigma is not small enough
    if f_alpha_k > 1e6
        f_alpha_k = 0;
        f_alpha_k1 = 0;
        f_alpha_k2 = 0;
        q = 0;
        sigma = sigma/10;
        disp('Golden Search sigma is decreased')
        i = i + 1;
        if i == termination
            disp('Golden search could not be converged. Stopping search algorithm...')
            optStepSize = NaN;
            return 
        end
    end

end

alpha_low = alpha_k;
alpha_up = alpha_k2;

%% phase 2 loop
I = alpha_up - alpha_low;
alpha_b = alpha_low + GR*I;
alpha_a = alpha_low + (1-GR)*I;

while I > eps

    f_alpha_b = subs(f_stepS,a,alpha_b);
    f_alpha_a = subs(f_stepS,a,alpha_a);

    if f_alpha_b > f_alpha_a
        alpha_up = alpha_b;
        alpha_b = alpha_a;
        I = alpha_up - alpha_low;
        alpha_a = alpha_low + (1-GR)*I;

    elseif f_alpha_b < f_alpha_a
        alpha_low = alpha_a;
        alpha_a = alpha_b;
        I = alpha_up - alpha_low;
        alpha_b = alpha_low + GR*I;
    
    else 
        alpha_low = alpha_a;
        alpha_up = alpha_b;
        I = alpha_up - alpha_low;
        alpha_b = alpha_low + GR*I;
        alpha_a = alpha_low + (1-GR)*I;
    end
end

optStepSize = (alpha_up + alpha_low) / 2;

end








