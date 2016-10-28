function Freproc = bockrisSwinkelTerm(theta,n)
    Freproc = theta./((1-theta).^n) .* (theta + n .* (1-theta)).^(n-1) ./ n.^n;
end