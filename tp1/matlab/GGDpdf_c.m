function px = GGDpdf_c(x,sig,bet,lam)

num= bet*(exp(-(x./sig).^bet)).*(x.^(bet*lam-1));
denom=gamma(lam)*sig^(bet*lam);

px = num/denom;