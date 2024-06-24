function [t,phid]=EulerCoupledReducedBeforeAverage(jj,samplenum,omega,phi0,h,ts,te,mu,outputTheandYandZy,JumpLtheta,JumpRtheta,JumpLthetaR,JumpRthetaL)

    N=floor((te-ts)/h);
    phi(1,:)=phi0;
    t(1)=0;
    yy=zeros(1,2);
    zy=zeros(1,2);
    indlr=zeros(1,2);

    for kk=1:2
        [~,lab]=min(abs(outputTheandYandZy(:,1)-phi(1,kk)));
        if outputTheandYandZy(lab,2)>=0
            indlr(kk)=2;
            if phi(1,kk)>=JumpRtheta(kk,1)
                phi(1,kk)=JumpRthetaL(kk,1);
                [~,lab]=min(abs(outputTheandYandZy(:,1)-phi(1,kk)));
                indlr(kk)=1;
            end    
        else
            indlr(kk)=1;
            if phi(1,kk)<=JumpLtheta(kk,1)
                phi(1,kk)=JumpLthetaR(kk,1);
                [~,lab]=min(abs(outputTheandYandZy(:,1)-phi(1,kk)));
                indlr(kk)=2;
            end   
        end
        yy(kk)=outputTheandYandZy(lab,3);
        zy(kk)=outputTheandYandZy(lab,4);
    end
    jjL=[jj*samplenum+1,jj*samplenum+1];
    jjR=[jj*samplenum+1,jj*samplenum+1];
    for ii=1:N
        t(ii+1)=t(ii)+h;
        phi(ii+1,1)=phi(ii,1)+(omega+mu*zy(1)*(yy(2)-yy(1)))*h;
        phi(ii+1,2)=phi(ii,2)+(omega+mu*zy(2)*(yy(1)-yy(2)))*h;
        for kk=1:2
            if indlr(kk)==1
                if (phi(ii+1,kk)>=JumpLtheta(kk,jjL(kk)+1))&&(phi(ii,kk)<JumpLtheta(kk,jjL(kk)+1))
                    phi(ii+1,kk)=JumpLthetaR(kk,jjL(kk)+1);
                    jjL(kk)=jjL(kk)+1;
                    indlr(kk)=2;
                end
            elseif indlr(kk)==2
                if (phi(ii+1,kk)>=JumpRtheta(kk,jjR(kk)+1))&&(phi(ii,kk)<JumpRtheta(kk,jjR(kk)+1))
                    phi(ii+1,kk)=JumpRthetaL(kk,jjR(kk)+1);
                    jjR(kk)=jjR(kk)+1;
                    indlr(kk)=1;
                end
            end
            phi(ii+1,kk)=mod(phi(ii+1,kk),2*pi);
            [~,lab]=min(abs(outputTheandYandZy(:,1)-phi(ii+1,kk)));
            yy(kk)=outputTheandYandZy(lab,3);
            zy(kk)=outputTheandYandZy(lab,4);
        end
    end
    phid=phi(:);
end
