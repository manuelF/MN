function [ imggg, newlabels ] = filterimages( imgs, labels, arein )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here    
    temp=zeros(size(imgs));
    templ=zeros(size(labels));
    thisindex=1;
    rows=size(imgs,1);
    for im=1:rows,
        
        if(any(ismember(arein,labels(im))))
            temp(thisindex,:)=imgs(im,:);
            templ(thisindex,:)=labels(im);
            thisindex=thisindex+1;
        end
    end
    imggg=temp(1:thisindex,:);
    newlabels=templ(1:thisindex,:);

end

