function [s1,s2]=delayshift(s1,s2,maxlag)
[S, l]=synch_DS(s1,s2,maxlag);
x1=1:length(s1);
x2=(1-l(1)):(length(s1)-l(1));
figure,plot(x1,s1),hold all,plot(x2,s2)