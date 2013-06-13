limit=10000;
width=784;

imgg=leerMNISTimage('Training Images',limit);
imgg=double(imgg);
X=imgg/sqrt(size(imgg,ndims(imgg))-1);
[U,S,V]=svd(X);
mu=mean(X,2);
M=zeros(limit,width);

for im=1:limit,
    M(im,:)=(X(im)-mu)'/sqrt(width-1);
end

[U,S,V]=svd(M);

xp=zeros(width,limit);
for im=1:limit,
xp(:,im)=V'*X(:,im);
end
scatter3(xp(:,1),xp(:,2),xp(:,3))
