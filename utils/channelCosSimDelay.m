function [cosSim]=channelCosSimDelay(V,lag)
[h,w]=size(V);
cosSim=zeros(h,h);

for i=1:h
    for j=1:h
        l=lag(i,j);
        
        cosSim(i,j)=getCosineSimilarity(V(i,:),circshift(V(j,:),l));
        
        if(isnan(cosSim(i,j)))
            
            disp('NAN')
        end
        
    end
end