limit=10000; %cuantas imagenes maximo leemos
width=784;

imggg=double(leerMNISTimage('Training Images',limit)); 
labels=leerMNISTlabel('Training Labels',limit);
hold on
for nrs=1:9,

    [imgg, nl] =filterimages(imggg,labels,[1]*nrs); %deja solo las iamgenes con label []

    nrows = size(imgg,1); %cuantas imagenes quedaron

    X=imgg/sqrt(nrows-1); %dividimos por la sqrt(n-1)
    mu=mean(X); %obtenemos la imagen media
    M=zeros(nrows ,width); %creamos una matriz temporal

    for im=1:nrows ,
        M(im,:)=(X(im)-mu); %obtenemos la varianza de cada imagen
    end
    M=M/sqrt(width-1); %dividimos por n-1
    
    [U,S,V]=svd(M); %sacamos la svd de esa

    % V = avec(M*M') covarianza

    xp=zeros(nrows ,width); %otra matriz temporal
    for im=1:nrows , 
        xp(im,:)=V'*X(im,:)'; %aplicamos la transformacion a cada elemento
    end
    scatter3(xp(:,1),xp(:,2),xp(:,3),3,nl) %ploteamos tres coordenadas
end

xpp=[xp,nl];