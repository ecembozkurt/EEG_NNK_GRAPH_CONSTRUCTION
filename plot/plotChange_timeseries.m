function [sumWeight, n_connect]= plotChange_timeseries(W)
[h,w]=size(W);
ch=h;
n_connect=zeros(ch,1);
sumWeight=zeros(ch,1);

for i=1:ch
    n_connect(i)=length(nonzeros(W(i,:)));%number of nonzero connections
    sumWeight(i)=sum(nonzeros(W(i,:)));
end

