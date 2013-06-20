
% Haciendo el TP
%%
%Leemos el input de imagenes y labels
from=1;
limit=1000; %cuantas imagenes maximo leemos
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

originalpoints=1000;

imgg=double(leerMNISTimage('Training Images',from,originalpoints)); 
labels=leerMNISTlabel('Training Labels',from,originalpoints);
%Obtenemos la transformada de cada imagen que usamos para la distancia
%de mas matchs
timgg=zeros(originalpoints,width);
for ti=1:originalpoints,
    timgg(ti,:)=(V*imgg(ti,:)')';
end
%%
%Corremos los casos de test
from=40000; %Desde esta imagen 
limit=40500; %hasta esta

test_imgs=double(leerMNISTimage('Training Images',from,limit)); 
test_labels=leerMNISTlabel('Training Labels',from,limit);
nrows=size(test_imgs,1);
[qty,n]=hist(double(test_labels),10); %qty por digito
upperk= 100; %cantidad de columnas como maximo que tomamos

alldistanceprogression=zeros(upperk,1); %Progreso distancia a todas (top100)
partialprogression=zeros(upperk,3); %Progreso por norma
progression=zeros(upperk,1); %Progreso bo3
metodonorma=[1,2,inf]; %las normas que usamos

alldistancedigitprogression=zeros(upperk,10); %Progreso por digito por distancia a todos

partialdigitprogression=zeros(upperk,3,10); %Progreso por norma por digito
digitprogression=zeros(upperk,10); %Progreso bo3 por digito

dists=zeros(originalpoints,2); % distancia a todos los puntos

alldists=1;

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
                alldistancedigitprogression(k,index)=1+alldistancedigitprogression(k,index);
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
            digitprogression(k,mejor_indice+1)=1+digitprogression(k,mejor_indice+1);
        end
        
    end

    %obtenemos los % para los metodos de reconocimento
    progression(k)=100*(hit/(nrows));
    partialprogression(k,:)=100*(partialhit(:)/(nrows));
    alldistanceprogression(k)=100*(hitdistance/(nrows));
    for nrm=1:3,        
            partialdigitprogression(k,nrm,:)=100*(partialdigitprogression(k,nrm,:)./qty(:));        
    end
    alldistancedigitprogression(k,:)=100*(alldistancedigitprogression(k,:)/qty(:));
    digitprogression(k,:)=100*(digitprogression(k,:)./qty(:));
end

%Vemos como evoluciona el hitrate a medida que k-> inf, comparamos
%data1 = usando solo norma2
%data2 = usando el mas repetido de las 3 normas
%data3 = usando la distancia a todos

%plot([partialprogression(:,2),progression(:),alldistanceprogression(:)]) 

%legend('Norma2','Mejor de 3 normas','Distancia Total','Location','SouthEast');
%xlabel('k cantidad de columnas');
%ylabel('% Hitrate');

normToUse=3;
%plot([partialdigitprogression(:,normToUse,1),partialdigitprogression(:,normToUse,2),partialdigitprogression(:,normToUse,3),partialdigitprogression(:,normToUse,4),partialdigitprogression(:,normToUse,5),partialdigitprogression(:,normToUse,6),partialdigitprogression(:,normToUse,7),partialdigitprogression(:,normToUse,8),partialdigitprogression(:,normToUse,9),partialdigitprogression(:,normToUse,10)])
%plot([digitprogression(:,1),digitprogression(:,2),digitprogression(:,3),digitprogression(:,4),digitprogression(:,5),digitprogression(:,6),digitprogression(:,7),digitprogression(:,8),digitprogression(:,9),digitprogression(:,10)])
plot([alldistancedigitprogression(:,1),alldistancedigitprogression(:,2),alldistancedigitprogression(:,3),alldistancedigitprogression(:,4),alldistancedigitprogression(:,5),alldistancedigitprogression(:,6),alldistancedigitprogression(:,7),alldistancedigitprogression(:,8),alldistancedigitprogression(:,9),alldistancedigitprogression(:,10)])
legend('[%] 0','[%] 1','[%] 2','[%] 3','[%] 4','[%] 5','[%] 6','[%] 7','[%] 8','[%] 9','Location','SouthEast');
xlabel('k cantidad de columnas');
ylabel('% Hitrate por digito');
title('Hitrate por digito de la deteccion por norma Infinito, para M Covarianza de 1k');


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
