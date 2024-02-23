function DTWDis=getDTW(V)
[h,w]=size(V);
for i=1:h
    for j=1:h
  %  V(i,:)=(V(i,:)-min(V(i,:)))/(max(V(i,:))-min(V(i,:)));
   % V(j,:)=(V(j,:)-min(V(j,:)))/(max(V(j,:))-min(V(j,:)));

    
   % V(i,:)=V(i,:)/norm(V(i,:));
   % V(j,:)=V(j,:)/norm(V(j,:));
    [DTWDis(i,j),~,~]=dtw(V(i,:),V(j,:));
    DTWDis(i,j)=1/(1+DTWDis(i,j));  %normalize
        
    end
end