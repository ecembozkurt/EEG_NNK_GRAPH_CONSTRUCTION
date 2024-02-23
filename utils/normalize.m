function [nA]=normalize(A)
nA=(A-min(A))/(max(A)-min(A));
