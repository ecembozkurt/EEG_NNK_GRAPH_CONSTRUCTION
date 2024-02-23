function [mutInf]=getMI(V)
[h,w]=size(V);
mutInf=zeros(h,h);
for i=1:h
    for j=1:h
        mutInf(i,j)=mi(V(i,:)',V(j,:)'); 
    end
end