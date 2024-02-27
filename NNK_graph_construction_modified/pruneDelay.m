function W=pruneDelay(W,maxDelay)
nolayer=maxDelay+1;
[CH,T]=size(W);
ch=CH/nolayer;
t=T/nolayer;
mask=zeros(CH,T);
mask2=zeros(CH,T);
mask3=ones(CH,T);

TTTT=W((1:ch),:)+W(((1+ch):2*ch),:)+W(((1+2*ch):3*ch),:);

for ich=1:ch
    for iii=1:T
        mask(ich+ch.*(-1+find( max(W([ich,ch+ich,2*ch+ich],iii),[],1))),iii)=1;
        mask3(ich,[ich,ch+ich,2*ch+ich])=0;
        mask3([ich,ch+ich,2*ch+ich],ich)=0;
    end
end

for ich=1:ch
    for iii=1:CH
        mask2(iii,ich+ch.*(-1+find( max(W(iii,[ich,ch+ich,2*ch+ich])))))=1;
    end
end
%figure,imagesc(mask)
%figure,imagesc(mask2)

TTTT=TTTT.*mask2(1:ch,:);
W1=W.*mask;
W2=W.*mask2;
W=W.*mask3;
mask_combined=min(mask,mask2);
%figure,imagesc(W1)
%figure,imagesc(W2)
W=W.*mask_combined;
W=W;