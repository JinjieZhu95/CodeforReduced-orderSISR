function [t,phi]=EulerPhi(phi0,h,ts,te,Gammathetad,delta,mu)

N=floor((te-ts)/h);
t(1)=ts;
phi(1)=phi0;
dd=2*pi/length(Gammathetad);
for ii=1:N
   t(ii+1)=t(ii)+h;
   
   phi0=mod(phi0,2*pi);
   I=floor(phi0/dd)+1;
   phi(ii+1)=phi(ii)+(delta+mu*Gammathetad(I))*h;
   phi0=phi(ii+1);
%    phi(ii+1)=mod(phi(ii+1),2*pi);

end