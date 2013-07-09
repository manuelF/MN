function [ V ] = calc_qr( M, its )
%CALCQR Calculate its iterations of qr eig method

    n=size(M,1);
    auVec=eye(n);
    for i=1:its,
        [Q,R]=qr(M);
        M=R*Q;
        auVec=auVec*Q;
    end

    auVal=diag(M);
    [~,indexes]=sort(auVal,'descend');
    V=zeros(size(auVec));
    for j=1:n,
        V(:,j)=auVec(:,indexes(j));
    end
end

