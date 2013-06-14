
% Haciendo el TP
from=1;
limit=40000; %cuantas imagenes maximo leemos
width=784;

imgg=double(leerMNISTimage('Training Images',from,limit)); 
labels=leerMNISTlabel('Training Labels',from,limit);

if ~(exist('covarianza.mat','file')==2)

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

allmeans=zeros(10,width);

for nro=0:9,
    [thisimgg,thislabels]=filterimages(imgg,labels,[1]*nro);
    nrows=size(thisimgg,1);
    M2=zeros(nrows,width);
    for im=1:nrows ,
        M2(im,:)=(V*thisimgg(im,:)')';
    end
    allmeans((nro+1),:)=mean(M2);
    
end

% ahora a probar esto

limit=40100;
from=40000;
test_imgs=double(leerMNISTimage('Training Images',from,limit)); 
test_labels=leerMNISTlabel('Training Labels',from,limit);
nrows=size(test_imgs);

k= 25; %cantidad de columnas que tomamos
for k=1:40,
    hit=0;
    for im=1:nrows,
        transf=(V*test_imgs(im,:)')'; % transformamos una imagen
        norma=1e16;
        mejor_indice=11;
        for comparar=1:10, %comparamos contra todos los mu-digitos
            q=norm(allmeans(comparar,1:k)-transf(1:k),2); %usando norma euclidia
            if(q<norma) %nos quedamos con la menor distancia
                norma=q;
                mejor_indice=comparar-1;
            end
        end
        str=sprintf('para el indice %d, se estimo %d y en realidad era %d\n',im,mejor_indice,test_labels(im));
        if(mejor_indice==test_labels(im)) %contabilizamos el hit
            hit=hit+1;
        end
        %disp(str);
    end
    str=sprintf('hitrate k=%d  : %.2f%%',k,100*(hit/(limit-from)));
    disp(str);
end


%xp=zeros(nrows ,width); %otra matriz temporal
%for im=1:nrows , 
%    xp(im,:)=V'*X(im,:)'; %aplicamos la transformacion a cada elemento
%end
%pointy=mean(xp);
%allmeans(nrs+1,:)=pointy(:);
%%scatter3(pointy(:,1),pointy(:,2),pointy(:,3),15,nl(1));
%scatter3(xp(:,1),xp(:,2),xp(:,3),3,nl) %ploteamos tres coordenadas

