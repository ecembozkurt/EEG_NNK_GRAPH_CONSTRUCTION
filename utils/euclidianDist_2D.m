function [D]=euclidianDist_2D(A,B)
D=sqrt(sum((A-B).^2));