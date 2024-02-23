function C=channelAVENEIGHDEG2d(VG)
ch_number=size(VG,2);
coef=zeros(ch_number,size(VG{1},1));
for i4=1:ch_number
    cd '../code'
    coef(i4,:) =ave_neighbor_deg(VG{i4});
    cd '../EEG-W-NNK-main'
end
C=coef;
