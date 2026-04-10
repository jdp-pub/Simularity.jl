normalized_legendre(s::Int) = [factorial(s)^2/factorial(2*s)*(-1)^(s-k)*binomial(s,k)*binomial(s+k,k) for k in 0:s]
polyroots(p::Vector) = eig_vals(hcat(I(length(p)-1)[:,2:length(p)-1],-p[1:length(p)-1]))

