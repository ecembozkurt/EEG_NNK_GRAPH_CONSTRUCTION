function N=interlayer_conn(W,ch_number)
N=0;
for i =1:size(W,1)
    for k=1:size(W,2)
        if floor((i-1)/ch_number)~=floor((k-1)/ch_number) && W(i,k)~=0
            % i~=k && 
            N= N+1;
           %rem((i-1),ch_number)+1,rem((k-1),ch_number)+1
        end
        
    end
end
N=N/2;