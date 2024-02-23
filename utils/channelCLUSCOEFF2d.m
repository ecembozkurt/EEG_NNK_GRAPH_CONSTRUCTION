function C=channelCLUSCOEFF2d(VG)
ch_number=size(VG,2);
coef=zeros(ch_number,size(VG{1},1));
for i4=1:ch_number
    [~, coef(i4,:)] = avgClusteringCoefficient(VG{i4});
end
C=coef;
