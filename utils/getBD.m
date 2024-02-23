function sim=getBD(P)
[h,w]=size(P);
sim=zeros(h,h);
for i=1:h
    for j=1:h
    sim(i,j)=BachDist(P(i,:),P(j,:));
    end
end