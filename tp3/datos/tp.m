limit=10000;
width=784;

imggg=double(leerMNISTimage('Training Images',limit)); 
labels=leerMNISTlabel('Training Labels',limit);


imgg=filterimages(imgg,labels,[0,3]);

nrows = size(imgg,1);

X=imgg/sqrt(size(imgg,ndims(imgg))-1);
%[U,S,V]=svd(X);
mu=mean(X);
M=zeros(nrows ,width);

for im=1:nrows ,
    M(im,:)=(X(im)-mu)/sqrt(width-1);
end

[U,S,V]=svd(M);

xp=zeros(nrows ,width);
for im=1:nrows ,
    xp(im,:)=V'*X(im,:)';
end
scatter3(xp(:,1),xp(:,2),xp(:,3))
