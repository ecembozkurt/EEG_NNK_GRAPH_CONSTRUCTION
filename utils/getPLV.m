function plv=getPLV(V)
[h,w]=size(V);
plv=zeros(h,h);
for i=1:h
    for j=1:h
    plv(i,j)=(1/w)*abs(sum( exp(1i*(phase(V(i,:)) -phase(V(j,:))))));     
    end
    
end
