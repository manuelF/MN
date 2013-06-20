
% Haciendo el TP
%%
%Leemos el input de imagenes y labels
from=1;
limit=30000; %cuantas imagenes maximo leemos
width=784;

imgg=double(leerMNISTimage('Training Images',from,limit)); 
labels=leerMNISTlabel('Training Labels',from,limit);
%%
%Calculamos una unica vez la matriz de covarianza

if ~(exist('covarianza.mat','file')==2)

    nrows = size(imgg,1); %cuantas imagenes quedaron

    mu=mean(imgg); %obtenemos la imagen media
    M=zeros(nrows ,width); %creamos una matriz temporal


    for im=1:nrows ,
        M(im,:)=(imgg(im)-mu); %obtenemos la varianza de cada imagen
    end
    M=M/sqrt(width-1); %dividimos por n-1

    [U,S,V]=svd(M); %sacamos la svd de esa
    clearvars U S M
    save('covarianza.mat','V');
end
load('covarianza.mat','V');
% V = avec(M*M') covarianza
%%
%Buscamos la media de todos las imagenes
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

originalpoints=200;
%Obtenemos la transformada de cada imagen que usamos para la distancia
%de mas matchs
timgg=zeros(originalpoints,width);
for ti=1:originalpoints,
    timgg(ti,:)=(V*imgg(ti,:)')';
end
%%
%Corremos los casos de test
from=00001; %Desde esta imagen 
limit=59000; %hasta esta

test_imgs=double(leerMNISTimage('Training Images',from,limit)); 
test_labels=leerMNISTlabel('Training Labels',from,limit);
nrows=size(test_imgs,1);
[qty,n]=hist(double(test_labels),10); %qty por digito
upperk= 100; %cantidad de columnas como maximo que tomamos

alldistanceprogression=zeros(upperk,1); %Progreso distancia a todas (top100)
partialprogression=zeros(upperk,3); %Progreso por norma
progression=zeros(upperk,1); %Progreso bo3
metodonorma=[1,2,inf]; %las normas que usamos

partialdigitprogression=zeros(upperk,3,10); %Progreso por norma por digito
digitprogression=zeros(upperk,1,10); %Progreso bo3 por digito

dists=zeros(originalpoints,2); % distancia a todos los puntos

alldists=0;

for k=1:upperk, %para cada k cantidad de columnas
    hit=0;
    hitdistance=0;
    partialhit=zeros(1,3);
    
    for im=1:nrows,
        transf=(V*test_imgs(im,:)')'; % transformamos una imagen
        if(alldists==1)
            %sacamos la distanacia norma2 a todas las imagenes de la V
            for id=1:originalpoints, 
                dists(id,1)=norm(timgg(id,1:k)-transf(1:k),2);            
            end
            %pongo las labels de los numeros apropiados
            dists(1:originalpoints,2)=labels(1:originalpoints);
            %ordenamos por menor distancia 
            dists=sortrows(dists);
            %nos quedamos con el top100
            dists=dists(1:200,:);
            %el histograma de los 10 digitos nos dice cuanto se repite c/u
            [appears]=hist(dists(:,2),10);
            %obtenemos el indice = digito que mas aparecio
            [elem,index]=max(appears);

            if((index-1)==test_labels(im)) %contabilizamos el hit si le acertamos
                hitdistance=hitdistance+1;
            end
        end
        
        %ahora hacemos las comparaciones de normas
        mejor_indice_p=[11,11,11];
               
        for inorma=1:3,
            norma=1e16;
            for comparar=1:10, %comparamos contra todos los mu-digitos
                q=norm(allmeans(comparar,1:k)-transf(1:k),metodonorma(inorma)); %usando norma euclidia
                if(q<norma) %nos quedamos con la menor distancia
                    norma=q;
                    mejor_indice_p(inorma)=comparar-1;
                end
            end
        end
        %nos quedamos ademas con el digito reconocido mas repetido, y defaultea en norma2
        
        mejor_indice=mejor_indice_p(2);
        if((mejor_indice_p(1)==mejor_indice_p(3))||(mejor_indice_p(1)==mejor_indice_p(2)))
            mejor_indice=mejor_indice_p(1);
        end
        if((mejor_indice_p(2)==mejor_indice_p(3))||(mejor_indice_p(1)==mejor_indice_p(2)))
            mejor_indice=mejor_indice_p(2);
        end
        %contabilizamos los hits parciales, si alguno le pego        
        for inorma=1:3,
            if (mejor_indice_p(inorma)==test_labels(im))
                partialhit(inorma)=partialhit(inorma)+1;
                partialdigitprogression(k,inorma,test_labels(im)+1)=1+partialdigitprogression(k,inorma,test_labels(im)+1);
            end
        end
        
        %contabilizamos el hit del bo3
        if(mejor_indice==test_labels(im)) 
            hit=hit+1;
        end
        
    end

    %obtenemos los % para los metodos de reconocimento
    progression(k)=100*(hit/(nrows));
    partialprogression(k,:)=100*(partialhit(:)/(nrows));
    alldistanceprogression(k)=100*(hitdistance/(nrows));
    for nrm=1:3,
        for dig=1:10,
            partialdigitprogression(k,nrm,dig)=100*(partialdigitprogression(k,nrm,dig)/qty(dig));
        end
    end
end

%Vemos como evoluciona el hitrate a medida que k-> inf, comparamos
%data1 = usando solo norma2
%data2 = usando el mas repetido de las 3 normas
%data3 = usando la distancia a todos

plot([partialprogression(:,2),progression(:),alldistanceprogression(:)]) 

legend('Norma2','Mejor de 3 normas','Distancia Total','Location','SouthEast');
xlabel('k cantidad de columnas');
ylabel('% Hitrate');


plot([partialdigitprogression(:,2,1),partialdigitprogression(:,2,2),partialdigitprogression(:,2,3),partialdigitprogression(:,2,4),partialdigitprogression(:,2,5),partialdigitprogression(:,2,6),partialdigitprogression(:,2,7),partialdigitprogression(:,2,8),partialdigitprogression(:,2,9),partialdigitprogression(:,2,10)])
legend('[%] 0','[%] 1','[%] 2','[%] 3','[%] 4','[%] 5','[%] 6','[%] 7','[%] 8','[%] 9','Location','SouthEast');
xlabel('k cantidad de columnas');
ylabel('% Hitrate por digito');
title('Hitrate de la deteccion por norma2');
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
