#!/usr/local/bin/octave -qf
format('long');;
args = argv();
baseOriginal = str2num(args{1});
longitudTotal = size(args,1)-2;
longitud=longitudTotal/2;
filename = args{end}
figure("visible","off");
originalLastNumber=baseOriginal+longitud;
x=baseOriginal:1:originalLastNumber-1;
#x=10.^(-x) #eleva todos a la 10^-x
disp(x);
cc=('brg');
vars=["Beta";"Iteraciones" ];
y=zeros(2,longitud)
 for q=0:1
    offset=q*longitud;
	lastnumber=baseOriginal+longitud+(offset);

	for i=1:longitud
		y(q+1,i)=str2num(args{offset+i+1});
	end
	disp(y);


#	hold on

#	plot(x, y,'-*','color',cc(q+1));

end
[haxes,hfig1, hfig2]=plotyy(x,[y(1,:)],x,[ y(2,:)]);
axes(haxes(1));
ylabel("Valor Estimado de Beta");
axes(haxes(2));
ylabel("Cantidad de iteraciones");
z=x;
z=10.^(-z);
set(haxes(1),'XTickLabel',[]);
set(haxes(2),'XTickLabel',z);
xlabel("Valor del epsilon");

set(hfig2,'LineStyle','--');

legend(vars,'location','NorthEast');
title(filename)
print(filename,"-dpng");
