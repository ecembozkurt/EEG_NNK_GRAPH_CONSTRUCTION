function plot_staticConn2(W_s,isknn,eeg_xyz,ch_number,st)
fig=figure;
p2=plot3(eeg_xyz(:,1),eeg_xyz(:,2),eeg_xyz(:,3),'r*') ;
for i=1:ch_number
    text(eeg_xyz(i,1),eeg_xyz(i,2),eeg_xyz(i,3)+1,st{i})
    hold all
    for k=1:ch_number
        if W_s(i,k)==0 || W_s(i,k)==Inf
        else
            cc= abs(W_s(i,k));
            p2=plot3([eeg_xyz(i,1),eeg_xyz(k,1)],[eeg_xyz(i,2),eeg_xyz(k,2)],[eeg_xyz(i,3),eeg_xyz(k,3)],'blue');
            grid on
            p2.LineWidth=cc*5;
        end
    end
end
hold off;
if isknn==0
    title('NNK')
else
    title('NNK')
end