function imagest2movies()
video = VideoWriter('newGfig_64.mp4','MPEG-4'); %create the video object
open(video); %open the file for writing
Files=dir('./Images2/noGfig_tw_*.*');
%video.FrameRate=1;
for k=1:length(Files)
 FileName=Files(k).name;
 filedir=sprintf('./Images2/%s',FileName);
 I = imread(filedir); %read the next image
 frame=im2frame(I);
 writeVideo(video,frame); %write the image to file
end
close(video); %close the file