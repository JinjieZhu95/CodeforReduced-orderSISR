function [t,phia]=par4originalphi(ss,xy1,xy2,h,ts,te,Sig,mu,xynDist,outputPositiveX,outputNegativeX)

[t,xy]=EulerCoupledOriginal(xy1,xy2,h,ts,te,Sig,mu);
xyN=xy(1:xynDist:end,:);

for ii=1:size(xyN,1)
    for kk=1:2
        if xyN(ii,2*kk-1)>=0
            [~,lab]=min(abs(xyN(ii,2*kk)-outputPositiveX(:,2)));
            phi(ii,kk)=outputPositiveX(lab,1);
        elseif xyN(ii,2*kk-1)<0
            [~,lab]=min(abs(xyN(ii,2*kk)-outputNegativeX(:,2)));
            phi(ii,kk)=outputNegativeX(lab,1);
        end
    end
end
phia=phi(:);