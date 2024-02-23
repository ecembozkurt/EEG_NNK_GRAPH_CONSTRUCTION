function [corrs]=getCorrelation(V,type)
[h,w]=size(V);
corrs=zeros(h,h);
for i=1:h
    for j=1:h
        
        if strcmp(type,'Pearson')
            corrs(i,j)=corr(V(i,:)',V(j,:)','Type','Pearson');
        elseif strcmp(type,'Spearman')
            corrs(i,j)=corr(V(i,:)',V(j,:)','Type','Spearman');
        elseif strcmp(type,'Kendall')
            corrs(i,j)=corr(V(i,:)',V(j,:)','Type','Kendall');
        else
           corrs(i,j)=corr(V(i,:)',V(j,:)'); 
        end
    end
end