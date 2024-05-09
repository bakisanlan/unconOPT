function alpha = sumalpha(n,sigma)

GR = (sqrt(5) - 1)/2;
alpha = 0;

for i=0:n
    alpha = alpha + sigma*(1/GR)^i;
end