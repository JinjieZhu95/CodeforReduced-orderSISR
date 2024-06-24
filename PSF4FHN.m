close,clc,clear

a=0.9;epsilon=0.0001;

x0=[0.1,0.1]; % initial condition for FHN
zxy=[1.0,1.0]; % initial condition for phase sensitivity

dt=1e-2*epsilon; % interation time for Euler method
Trel=200; % relaxation time for limit cycle
Rep=2; % relaxtion time for adjoint equation
Tstep=1; % time step for output, i.e., Tstep*dt;
Dx=1e-6; Dy=1e-6; % stepsize for x,y derivatives

% Discard initial relaxation to the limit cycle orbit
for i=1:Trel
    X=evolX(x0,dt,epsilon,a);
    x0=X;
end

% Look for the phase origin,where first variable crosses zero from below
x1=x0(2);
while 1
    ox=x1;
    x0=evolX(x0,dt,epsilon,a);
    x1=x0(2);
    if ox<0 && x1>0
        break;
    end
end

ttp=0; % measure the oscillation period
while 1
    x1=x0(2);
    x0=evolX(x0,dt,epsilon,a);
    ttp=ttp+1;
    if x1<0 && x0(2)>0
        break;
    end
end

T=ttp*dt;omega=2*pi/T;

% record one period
xyrec=zeros(ttp,2);
xyrec(1,:)=x0;
for i=1:ttp-1
    xyrec(i+1,:)=evolX(xyrec(i,:),dt,epsilon,a);
end

% integrate the adjoint equation
% % output=zeros(1,5); % theta,x,y,zx,zy
% % num=0;
% % for tt=ttp:-1:1
% %     for i=1:Rep
% %         x0=xyrec(tt-i+1,:);
% %         zxy=evolZ(zxy,x0,I,epsilon,c,d,omega,Dx,Dy,dt);
% %     end
% %     if mod(tt,Tstep)==0
% %         num=num+1;
% %         output(num,:)=[omega*(tt-Rep+1)*dt,x0,zxy];
% %     end
% % end

output=zeros(1,5); % theta,x,y,zx,zy
num=0;
for i=1:Rep  %%%% after test, Rep=2 is enough for convergence
    for tt=ttp:-1:1
        x0=xyrec(tt,:);
        zxy=evolZ(zxy,x0,epsilon,a,omega,Dx,Dy,dt);
        
        if (i==Rep)&&(mod(tt,Tstep)==0)
        num=num+1;
        output(num,:)=[omega*tt*dt,x0,zxy];
        end
    end
end

% set(gca,'XTick',[0:pi:2*pi])
% set(gca,'xtickLabel',{'0','¦Ð','2¦Ð'})
        
    