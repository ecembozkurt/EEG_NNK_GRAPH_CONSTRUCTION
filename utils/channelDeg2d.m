function D=channelDeg2d(A)
ch_number=size(A,2);
T=size(A{1},1);
D=zeros(ch_number,T);
for ii =1:T
    for i=1:ch_number
        D(i,ii)=length(nonzeros(A{i}(ii,:)));
    end
end