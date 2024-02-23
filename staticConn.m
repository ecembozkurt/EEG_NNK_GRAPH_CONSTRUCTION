function W_s=staticConn(W)
n_wind= size(W,1);
ch=size(W,2);
ch_d=size(W,3);
W_s=ones(ch,ch_d);
for i= 1:(n_wind-1)
  W_s=squeeze(W(i,:,:)>0) & W_s;
end
