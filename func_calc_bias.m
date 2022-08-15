function [bs, MAE, logbs, logMAE]=func_calc_bias(M_in,O_in)

idx= isfinite(M_in) & isfinite(O_in);
M=M_in(idx);
O=O_in(idx);
N=length(M);

bs=(1./N).*sum(M-O);
MAE=(1./N).*sum(abs(M-O));

idx= isfinite(log10(M_in)) & isfinite(log10(O_in));
M=M_in(idx);
O=O_in(idx);
N=length(M);

logbs=10.^((1./N).*sum(log10(M)-log10(O)));
logMAE=10.^((1./N).*sum(abs(log10(M)-log10(O))));


end