figure, subplot 311, plot(annotat(1,:),'LineWidth',2), title('Expert 1')
subplot 312, plot(annotat(2,:), 'LineWidth',2), title('Expert 2')
subplot 313, plot(annotat(3,:),'LineWidth',2) ,title('Expert 3')
sgtitle('Annotations by Experts')

inrped_annot3_5=zeros(3,length(EEG5(1,:)));
inrped_annot3_5(1,:)=interp1(annotat5(1,:),(1/256):(1/256):length(annotat5(1,:)),'linear');
inrped_annot3_5(2,:)=interp1(annotat5(2,:),(1/256):(1/256):length(annotat5(1,:)),'linear');
inrped_annot3_5(3,:)=interp1(annotat5(3,:),(1/256):(1/256):length(annotat5(1,:)),'linear');

figure, subplot 311, plot(inrped_annot3_23(1,:),'LineWidth',2), title('Expert 1')
subplot 312, plot(inrped_annot3_23(2,:), 'LineWidth',2), title('Expert 2')
subplot 313, plot(inrped_annot3_23(3,:),'LineWidth',2) ,title('Expert 3')
sgtitle('Annotations by Experts')