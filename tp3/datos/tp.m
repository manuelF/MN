
% Haciendo el TP
from=1;
limit=10000; %cuantas imagenes maximo leemos
width=784;

imgg=double(leerMNISTimage('Training Images',from,limit)); 
labels=leerMNISTlabel('Training Labels',from,limit);

%calculamos una unica vez la matriz de covarianza

if ~(exist('covarianza.mat','file')==2)

    nrows = size(imgg,1); %cuantas imagenes quedaron

    mu=mean(X); %obtenemos la imagen media
    M=zeros(nrows ,width); %creamos una matriz temporal


    for im=1:nrows ,
        M(im,:)=(X(im)-mu); %obtenemos la varianza de cada imagen
    end
    M=M/sqrt(width-1); %dividimos por n-1

    [U,S,V]=svd(M); %sacamos la svd de esa
    %clearvars U S M
    save('covarianza.mat','V');
end
load('covarianza.mat','V');
% V = avec(M*M') covarianza

allmeans=zeros(10,width); % vector que va a tener el promedio de los digs transformados

for nro=0:9,
    %obtenemos solo las imagenes que corresponden a este digito
    [thisimgg,thislabels]=filterimages(imgg,labels,[1]*nro); 
    nrows=size(thisimgg,1);
    M2=zeros(nrows,width);
    %transformamos todas
    for im=1:nrows ,
        M2(im,:)=(V*thisimgg(im,:)')';
    end
    %guardamos la imagen promedio
    allmeans((nro+1),:)=mean(M2);
    
end

% ahora a probar esto

limit=40100;
from=40000;
test_imgs=double(leerMNISTimage('Training Images',from,limit)); 
test_labels=leerMNISTlabel('Training Labels',from,limit);
nrows=size(test_imgs);

upperk= width; %cantidad de columnas que tomamos
progression=zeros(upperk); %anotamos para ver el progreso
for k=1:upperk,
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
    %disp(str);
    progression(k)=100*(hit/(limit-from));
end

plot(progression); %vemos como evoluciono el hitrate k-> inf

%Formas de medir:
%
%   Tomamos los t (parametro, o fijo: 100) mas cercanos.
%   Elegimos el digito mas repetido.
%   Si hay mas de uno tomamos los 2t.
%   3t, 4t... asi hasta que no haya repetidos.
%   Tambien podemos iterar 2t,3t,4t.. hasta que no haya repetidos y un
%   digito aparezca al menos el 20% de las veces (no se si esto sera muy
%   bueno).
%   Para medir distancia medimos distintas normas (norma 1, norma 2, norma
%   infinito).
%   
%
%
