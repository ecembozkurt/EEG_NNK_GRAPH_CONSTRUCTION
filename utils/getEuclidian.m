function EucDis=getEuclidian(V)
[h,w]=size(V);
EucDis=zeros(h,h);

for i=1:h
       %V(i,:)= V(i,:)-mean(V(i,:));
    
end
for i=1:h
    for j=1:h
  %  V(i,:)=(V(i,:)-min(V(i,:)))/(max(V(i,:))-min(V(i,:)));
   % V(j,:)=(V(j,:)-min(V(j,:)))/(max(V(j,:))-min(V(j,:)));

    
   % V(i,:)=V(i,:)/norm(V(i,:));
   % V(j,:)=V(j,:)/norm(V(j,:));
    EucDis(i,j)=euclidianDist_2D(V(i,:),V(j,:));
        
    end
end