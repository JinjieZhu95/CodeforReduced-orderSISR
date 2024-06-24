function [t,xy]=EulerCoupledOriginal(xy1,xy2,h,ts,te,Sig,mu)
a=0.9;epsilon=0.0001;
    N=floor((te-ts)/h);
    xy=zeros(N,4);
    xy(1,:)=[xy1,xy2];
    t(1)=0;
    for ii=1:N
        t(ii+1)=t(ii)+h;
        xy(ii+1,1)=xy(ii,1)+(xy(ii,1)-xy(ii,1).^3/3-xy(ii,2))/epsilon*h+sqrt(h*Sig(1)/epsilon)*randn(1);
        xy(ii+1,2)=xy(ii,2)+h*(xy(ii,1)+a)+mu*(xy(ii,4)-xy(ii,2))*h;
        xy(ii+1,3)=xy(ii,3)+(xy(ii,3)-xy(ii,3).^3/3-xy(ii,4))/epsilon*h+sqrt(h*Sig(2)/epsilon)*randn(1);
        xy(ii+1,4)=xy(ii,4)+h*(xy(ii,3)+a)+mu*(xy(ii,2)-xy(ii,4))*h;
    end
end