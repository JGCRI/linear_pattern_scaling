function varout=nanfill2d(varin)

varout=varin;
I=size(varin);
    Q=sum(isnan(varin)==1,1);
    QI=find(Q>=I(1)/2);
    L=length(QI);
    if L>0;
        L1=QI(1:L/2);
        L2=QI(L/2+1:end);
        for j=1:L/2;
            varout(:,j)=varout(:,L1(end)+1);
            varout(:,I(2)+1-j)=varout(:,L2(1)-1);
        end
    end  
    R=sum(isnan(varout)==1,2);
    RI=find(R>=I(2)/2);
    M=length(RI);
    I2=RI(RI>I(1)/2);
    I1=RI(RI<I(1)/2);
    if length(I2)>0;
        lastfull=varout(I2(1)-1,:);
    else
        lastfull=varout(end,:);
    end
    if length(I1)>0;
        firstfull=varout(I1(end)+1,:);
    else
        firstfull=varout(1,:);
    end
    if length(I1)>0;
        for j=1:length(I1);
            varout(j,:)=firstfull*(0.5+(j-1)/length(I1))+lastfull*(0.5-(j-1)/length(I1));
        end
    end
    if length(I2)>0;
        for j=1:length(I2);
            varout(I(1)+1-j,:)=lastfull*(0.5+(j-1)/length(I2))+firstfull*(0.5-(j-1)/length(I2));
        end
    end
end