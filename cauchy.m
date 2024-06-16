% Cauchy Mutation Strategy...
% This code was written bu U.Yuzgec on October 2022

function [output]=cauchy(X)
	n=size(X,2);
    output=zeros(size(X));
	for i=1:n
		output(:,i)=X(:,i)*(1+tan(pi*(rand-0.5))); 
	end
end
