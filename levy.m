% Levy Mutation Strategy...
% This code was written bu U.Yuzgec on October 2022

function [output]=levy(X)
	n=size(X,2);
	beta=3/2;
	sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    % output=zeros(size(X));
	u=randn(1,n)*sigma;
	v=randn(1,n);
	output_levy=u./abs(v).^(1/beta);
	output=X.*(1+output_levy);
end