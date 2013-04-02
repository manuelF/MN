#!/usr/local/bin/octave -qf
args = argv();
baseOriginal = str2num(args{1});
longitudTotal = size(args,1)-2;
longitud=longitudTotal/3;
filename = args{end}
figure("visible","off");
originalLastNumber=baseOriginal+longitud;
x=baseOriginal:1:originalLastNumber-1;
cc=lines(12);
vars=["Beta";"Lambda";"Sigma" ];
for q=0:2
	offset=q*longitud;
	disp(offset);
	lastnumber=baseOriginal+longitud+(offset);


	y=[1:1:longitud]

	for i=1:longitud
		y(i)=str2num(args{offset+i+1});
	end
	disp(y);


	hold on

	plot(x, y,'-*','color',cc((q+1),:));


end
ylabel("Valor de las variables");
xlabel("Digitos de precision");
legend(vars,'location','NorthEast');
title(filename)
print(filename,"-dpng");
