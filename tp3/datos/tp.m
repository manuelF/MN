
% Haciendo el TP
from=1;
limit=6000; %cuantas imagenes maximo leemos
width=784;

if ~(exist('covarianza.mat','file')==2)

    imgg=double(leerMNISTimage('Training Images',from,limit)); 
    labels=leerMNISTlabel('Training Labels',from,limit);
    allmeans=zeros(10,width);


    nrows = size(imgg,1); %cuantas imagenes quedaron

    X=imgg/sqrt(nrows-1); %dividimos por la sqrt(n-1)
    mu=mean(X); %obtenemos la imagen media
    M=zeros(nrows ,width); %creamos una matriz temporal


    for im=1:nrows ,
        M(im,:)=(X(im)-mu); %obtenemos la varianza de cada imagen
    end
    M=M/sqrt(width-1); %dividimos por n-1

    [U,S,V]=svd(M); %sacamos la svd de esa
    save('covarianza.mat','V');
end
load('covarianza.mat','V');
% V = avec(M*M') covarianza

xp=zeros(nrows ,width); %otra matriz temporal
for im=1:nrows , 
    xp(im,:)=V'*X(im,:)'; %aplicamos la transformacion a cada elemento
end
pointy=mean(xp);
allmeans(nrs+1,:)=pointy(:);
%scatter3(pointy(:,1),pointy(:,2),pointy(:,3),15,nl(1));
scatter3(xp(:,1),xp(:,2),xp(:,3),3,nl) %ploteamos tres coordenadas

