function e = Sh_Entropy(x)
p_x = x.^2./sum(x.^2);
p_x = p_x(find(p_x));
e = sum( -p_x.*log(p_x) );