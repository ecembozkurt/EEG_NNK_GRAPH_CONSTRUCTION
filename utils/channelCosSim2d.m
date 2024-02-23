function D=channelCosSim2d(A)
ch_number=size(A,2);
D=zeros(ch_number,ch_number);
for ii =1:ch_number
    for i=1:ch_number
        D(i,ii)=sum(sum(pdist2(A{i},A{ii},'cosine')));
    end
end