#!/usr/local/bin/octave -qf
args = argv();
base = str2num(args{1});
longitud = size(args,1)-3;
filename = args{end}
lastnumber=base+longitud;

x=base:1:lastnumber;
y=[1:1:longitud+1]

for i=1:longitud+1
	y(i)=str2num(args{i+1});
end

#disp(x);
class(y);
disp(y);


figure("visible","off");
hold on

plot(x, y,'-*');
title(filename)
print(filename,"-dpng");