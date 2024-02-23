function [wg]=GaussianKernel (D, p)
% d: distance
% p: reg. parameter
[h,w]=size(D);
wg=zeros(h,w);
for i=1:h
    for j=1: w
        wg(i,j)=exp(-D(i,j)^2/(2*p^2));
        
    end
end