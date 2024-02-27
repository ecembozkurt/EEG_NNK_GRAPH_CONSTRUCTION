cd '/Volumes/Seagate Backup Plus Drive/EEG SImulation/fieldtrip-20220607/plotting'
load vol.mat
figure(1450);
ft_plot_mesh(vol.bnd(1), 'edgecolor','none','facealpha',0.5,'facecolor',[0.6 0.6 0.8]);
hold all
plot3(eeg_xyz(:,2)+20,eeg_xyz(:,1),eeg_xyz(:,3)+50,'r*','LineWidth',2) ;
for i=1:ch_number
    text(eeg_xyz(i,2)+20,eeg_xyz(i,1),eeg_xyz(i,3)+51,st{i})
    hold all
    for k=1:ch_number
        if W_knn(i,k)==0 || W_knn(i,k)==Inf
        else
            cc= abs(W_knn(i,k));
            p2=plot3([eeg_xyz(i,2)+20,eeg_xyz(k,2)+20],[eeg_xyz(i,1),eeg_xyz(k,1)],[eeg_xyz(i,3)+50,eeg_xyz(k,3)+50],'blue');
            grid on
            p2.LineWidth=cc*5;
            
        end
    end
end
hold off;
title('KNN')


figure(1451);
ft_plot_mesh(vol.bnd(1), 'edgecolor','none','facealpha',0.5,'facecolor',[0.6 0.6 0.8]);
hold all
plot3(eeg_xyz(:,2)+20,eeg_xyz(:,1),eeg_xyz(:,3)+50,'r*','LineWidth',2) ;
for i=1:ch_number
    text(eeg_xyz(i,2)+20,eeg_xyz(i,1),eeg_xyz(i,3)+51,st{i})
    hold all
    for k=1:ch_number
        if W_nnk(i,k)==0 || W_nnk(i,k)==Inf
        else
            cc= abs(W_nnk(i,k));
            p2=plot3([eeg_xyz(i,2)+20,eeg_xyz(k,2)+20],[eeg_xyz(i,1),eeg_xyz(k,1)],[eeg_xyz(i,3)+50,eeg_xyz(k,3)+50],'blue');
            grid on
            p2.LineWidth=cc*5;
            
        end
    end
end
hold off;
title('NNK')
