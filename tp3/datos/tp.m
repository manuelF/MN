
% Haciendo el TP
%% Reading Input
%Leemos el input de imagenes y labels
from=1;
limit=15000; %cuantas imagenes maximo leemos
width=784;

imgg=double(leerMNISTimage('Training Images',from,limit)); 
labels=leerMNISTlabel('Training Labels',from,limit);
%% Covariance Calculation
%Calculamos una unica vez la matriz de covarianza

if ~(exist('covarianza.mat','file')==2)

    nrows = size(imgg,1); %cuantas imagenes quedaron

    mu=mean(imgg); %obtenemos la imagen media
    M=zeros(nrows ,width); %creamos una matriz temporal


    for im=1:nrows ,
        M(im,:)=(imgg(im)-mu); %obtenemos la varianza de cada imagen
    end
    M=M/sqrt(width-1); %dividimos por n-1
    M=M'*M;
    V=calc_qr(M,200);


   
    %[~,~,V]=svd(M); %sacamos la svd de esa
    %clearvars U S M
    %save('covarianza.mat','V');
end
%load('covarianza.mat','V');
V=V';
% V = avec(M*M') covarianza
%% Norm Average
%Buscamos la media de todos las imagenes
allmeans=zeros(10,width); % vector que va a tener el promedio de los digs transformados

for nro=0:9,
    %obtenemos solo las imagenes que corresponden a este digito
    [thisimgg,~]=filterimages(imgg,labels,[1]*nro); 
    nrows=size(thisimgg,1);
    M2=zeros(nrows,width);
    %transformamos todas
    for im=1:nrows ,
        M2(im,:)=(V*thisimgg(im,:)')';
    end
    %guardamos la imagen promedio
    allmeans((nro+1),:)=mean(M2);
    
end

%% Top 100 calculate and store 
originalpoints=1000;

imgg=double(leerMNISTimage('Training Images',from,originalpoints)); 
labels=leerMNISTlabel('Training Labels',from,originalpoints);
%Obtenemos la transformada de cada imagen que usamos para la distancia
%de mas matchs
timgg=zeros(originalpoints,width);
for ti=1:originalpoints,
    timgg(ti,:)=(V*imgg(ti,:)')';
end

%% Detection
%Corremos los casos de test
from=30000; %Desde esta imagen 
limit=30500; %hasta esta
test_imgs=double(leerMNISTimage('Training Images',from,limit)); 
test_labels=leerMNISTlabel('Training Labels',from,limit);
nrows=size(test_imgs,1);
[qty,~]=hist(double(test_labels),10); %qty por digito
upperk= 100; %cantidad de columnas como maximo que tomamos
fromk = 100;
alldistanceprogression=zeros(upperk,1); %Progreso distancia a todas (top100)
partialprogression=zeros(upperk,3); %Progreso por norma
progression=zeros(upperk,1); %Progreso bo3
metodonorma=[1,2,inf]; %las normas que usamos

alldistancedigitprogression=zeros(upperk,10); %Progreso por digito por distancia a todos

partialdigitprogression=zeros(upperk,3,10); %Progreso por norma por digito
digitprogression=zeros(upperk,10); %Progreso bo3 por digito

dists=zeros(originalpoints,2); % distancia a todos los puntos

detectedvalue=zeros(4,upperk,10,10); %Que valor detecté y que valor era correcto p/cada k

alldists=0;
normas=1;

for k=fromk:upperk, %para cada k cantidad de columnas
    hit=0;
    hitdistance=0;
    partialhit=zeros(1,3);
    
    for im=1:nrows,
        tmp=(V*test_imgs(im,:)')'; % transformamos una imagen        
        transf=tmp(1:k); %la recortamos a k columnas
        
        if(alldists==1)
            %sacamos la distanacia norma2 a todas las imagenes de la V
            for id=1:originalpoints, 
                dists(id,1)=norm(timgg(id,1:k)-transf,2);            
            end
            %pongo las labels de los numeros apropiados
            dists(1:originalpoints,2)=labels(1:originalpoints);
            %ordenamos por menor distancia 
            dists=sortrows(dists);
            %nos quedamos con el top200
            dists=dists(1:200,:);
            %el histograma de los 10 digitos nos dice cuanto se repite c/u
            [appears]=hist(dists(:,2),10);
            %obtenemos el indice = digito que mas aparecio
            [~,index]=max(appears);

            detectedvalue(4,k,index,test_labels(im)+1)=1+detectedvalue(4,k,index,test_labels(im)+1);
            if((index-1)==test_labels(im)) %contabilizamos el hit si le acertamos
                hitdistance=hitdistance+1;
                alldistancedigitprogression(k,index)=1+alldistancedigitprogression(k,index);
            end
        end
        
        %ahora hacemos las comparaciones de normas
        if(normas==1)
            mejor_indice_p=[11,11,11];

            mediaparcial=allmeans(:,1:k);

            for inorma=1:3,
                norma=1e16;
                for comparar=1:10, %comparamos contra todos los mu-digitos
                    q=norm(mediaparcial(comparar,:)-transf,metodonorma(inorma)); %usando norma euclidia
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
                detectedvalue(inorma,k,mejor_indice_p(inorma)+1,test_labels(im)+1)=1+detectedvalue(inorma,k,mejor_indice_p(inorma)+1,test_labels(im)+1);
            
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
        
    end

    %obtenemos los % para los metodos de reconocimento
    progression(k)=100*(hit/(nrows));
    partialprogression(k,:)=100*(partialhit(:)/(nrows));
    alldistanceprogression(k)=100*(hitdistance/(nrows));
    for  nrm=1:3,
        for dig=1:10,
            partialdigitprogression(k,nrm,dig)=100*(partialdigitprogression(k,nrm,dig)/qty(dig));        
        end
    end
    
    for dig=1:10,
        alldistancedigitprogression(k,dig)=100*(alldistancedigitprogression(k,dig)/qty(dig));
        digitprogression(k,dig)=100*(digitprogression(k,dig)/qty(dig));
    end
end
disp(partialprogression(100,2));
error('done');
%% Plotting
%Vemos como evoluciona el hitrate a medida que k-> inf, comparamos
%data1 = usando solo norma2
%data2 = usando el mas repetido de las 3 normas
%data3 = usando la distancia a todos
% 
h=figure;
 plot([partialprogression(:,2),progression(:),alldistanceprogression(:)]) 
% 
 legend('Norma2','Mejor de 3 normas','Top 100 Menor Dist','Location','SouthEast');
 xlabel('k cantidad de columnas','FontSize',14.0);
 ylabel('% Hitrate','FontSize',14.0);
 title('Hitrate por distintos metodos de deteccion, para M Covarianza de 30k','FontSize',14.0);
 
 
 
 
filename=sprintf('../informe/plots/hitrate-30kcv');
    saveas(h,filename,'png')
    saveas(h,filename,'fig')
 error('esto no');
for im=1:3,
    h=figure;
    normNombre={'norma1','norma2','normaInf','top100','BestOf3'};
    normToUse=im;
    %plot([partialdigitprogression(:,normToUse,1),partialdigitprogression(:,normToUse,2),partialdigitprogression(:,normToUse,3),partialdigitprogression(:,normToUse,4),partialdigitprogression(:,normToUse,5),partialdigitprogression(:,normToUse,6),partialdigitprogression(:,normToUse,7),partialdigitprogression(:,normToUse,8),partialdigitprogression(:,normToUse,9),partialdigitprogression(:,normToUse,10)])
    plot([digitprogression(:,1),digitprogression(:,2),digitprogression(:,3),digitprogression(:,4),digitprogression(:,5),digitprogression(:,6),digitprogression(:,7),digitprogression(:,8),digitprogression(:,9),digitprogression(:,10)])

    %plot([alldistancedigitprogression(:,1),alldistancedigitprogression(:,2),alldistancedigitprogression(:,3),alldistancedigitprogression(:,4),alldistancedigitprogression(:,5),alldistancedigitprogression(:,6),alldistancedigitprogression(:,7),alldistancedigitprogression(:,8),alldistancedigitprogression(:,9),alldistancedigitprogression(:,10)])
    legend('[%] 0','[%] 1','[%] 2','[%] 3','[%] 4','[%] 5','[%] 6','[%] 7','[%] 8','[%] 9','Location','SouthEast');
    xlabel('k cantidad de columnas','FontSize',14.0);
    ylabel('% Hitrate por digito','FontSize',14.0);
    title(sprintf('Hitrate por digito de la deteccion por %s, para M Covarianza de 30k',normNombre{normToUse}),'FontSize',14.0);
    filename=sprintf('../informe/plots/pordig-30kcv-%s',normNombre{normToUse});
    saveas(h,filename,'png')
    saveas(h,filename,'fig')
end


%%%%%%HEAT MAP DE EQUIVOCACIONES
%%%%% EJE X = detecté
%%%%% EJE Y = valor posta

digitos={'0','1','2','3','4','5','6','7','8','9'};
normNombre={'norma_1','norma_2','norma_inf','top100'};
toplot=[5,10,25,50,100];
for im=1:4,
    im=5;
    h=figure;
    kelegido=toplot(im);
    normToUse=2;
    colormap('gray');
    imagesc(-reshape(detectedvalue(normToUse,kelegido,:,:),[10 10]));
    set(gca,'dataAspectRatio',[1 1 1]);
    set(gca,'YTick',[1:10]);
    set(gca,'XTick',[1:10]);
    set(gca,'YTickLabel',digitos);
    set(gca,'XTickLabel',digitos);
    xlabel('Valor detectado');
    ylabel('Valor verdadero');
    titulo=sprintf('Cantidad de hits por digito con deteccion usando %s, con K=%d, M Covarianza 30k',normNombre{normToUse},kelegido);
    title(titulo);

    colorbar;

    filename=sprintf('../informe/plots/heatmap-30kcv-k%d-%s',kelegido,normNombre{normToUse});
    saveas(h,filename,'png')
    saveas(h,filename,'fig')
end