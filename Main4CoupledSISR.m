clear,clc,close

%% (1) PDFy MLE for sigma=0.01 and 0.012
load('LinearAndThirdOrderFitOfAandBversusSigma4Weibull.mat') %%obtained from the code file 'Estimateofabversussigma.m'
Sig=[0.01,0.012];

YL=-0.66:0.001:-0.3;
aL=pla(1)*Sig+pla(2);
bL=plb(1)*Sig.^3+plb(2)*Sig.^2+plb(3)*Sig+plb(4);

EstimateWLeft=zeros(length(Sig),length(YL));
for ii=1:length(Sig)
    EstimateWLeft(ii,:)=bL(ii)./aL(ii).*(abs(YL)./aL(ii)).^(bL(ii)-1).*exp(-(abs(YL)./aL(ii)).^bL(ii));
end

YR=0.3:0.001:0.67;
aR=pra(1)*Sig+pra(2);
bR=prb(1)*Sig.^3+prb(2)*Sig.^2+prb(3)*Sig+prb(4);
EstimateWRight=zeros(length(Sig),length(YR));
for ii=1:length(Sig)
    EstimateWRight(ii,:)=bR(ii)./aR(ii).*(abs(YR)./aR(ii)).^(bR(ii)-1).*exp(-(abs(YR)./aR(ii)).^bR(ii));
end
plot(YL,EstimateWLeft(1,:))
hold on
plot(YR,EstimateWRight(1,:))
%% (2) PDFtheta sigma=0.01 and 0.012
% load('Sigma0d1and0d12PDFy.mat') %%% obtained from the above (1)
% load('Zxya0d9epsilon1em4dt1em6.mat') %% Phase sensitivity function result obtained by the adjoint method. The adjoint method can be referred from e.g., Nakao 2016 Contemporary Physics and code therein.
% %%%%%%%%%%%% left %%%%%%%%%%%%%%%%%%%%%
% labelL=find(output(:,2)<0);
% theL=output(labelL,1);YtheL=output(labelL,3);
% labL=zeros(1,length(YL));
% for ii=1:length(YL)
%     [~,labL(ii)]=min(abs(YL(ii)-YtheL));
% end
% TheL=theL(labL);
% normleft1=abs(trapz(TheL,EstimateWLeft(1,:)));
% normleft2=abs(trapz(TheL,EstimateWLeft(2,:)));
% 
% %%%%%%%%%%%% right %%%%%%%%%%%%%%%%%%%%%
% labelR=find(output(:,2)>0);
% theR=output(labelR,1);YtheR=output(labelR,3);
% labR=zeros(1,length(YR));
% for ii=1:length(YR)
%     [~,labR(ii)]=min(abs(YR(ii)-YtheR));
% end
% TheR=theR(labR);
% normright1=abs(trapz(TheR,EstimateWRight(1,:)));
% normright2=abs(trapz(TheR,EstimateWRight(2,:)));

%% (3) Ktheta versus PDFtheta
% load('Ktheta(JumpStrength).mat') %% obtained from the code file 'ThetaJumpleft2right.m'
% load('Sigma0d1and0d12PDFtheta.mat') %%% obtained from the above (2)
% thelK=thetal(labL);
% therK=thetar(labR);
% LL=zeros(1,length(TheL));
% LR=zeros(1,length(TheR));
% for ii=1:length(TheL)
%     [~,LL(ii)]=min(abs(TheL(ii)-thelK));
% end
% for ii=1:length(TheR)
%     [~,LR(ii)]=min(abs(TheR(ii)-therK));
% end
% TheKL=KthetaL(LL);TheKR=KthetaR(LR);
% Kthe1=abs(trapz(TheL,TheKL'.*EstimateWLeft(1,:)/normleft1))+abs(trapz(TheR,TheKR'.*EstimateWRight(1,:)/normright1));
% Kthe2=abs(trapz(TheL,TheKL'.*EstimateWLeft(2,:)/normleft2))+abs(trapz(TheR,TheKR'.*EstimateWRight(2,:)/normright2));

%% (4) Coupling function Gammad(phi)
%%%%%%%%%%%%%%%%% mean theta on the left and right %%%%%%%%%%%%%%
% load('Sigma0d1and0d12PDFtheta.mat') %%Obtained from the above (2)
% load('LinearAndThirdOrderFitOfAandBversusSigma4Weibull.mat') %%obtained from the code file 'Estimateofabversussigma.m'
% aL=pla(1)*Sig+pla(2);
% bL=plb(1)*Sig.^3+plb(2)*Sig.^2+plb(3)*Sig+plb(4);
% aR=pra(1)*Sig+pra(2);
% bR=prb(1)*Sig.^3+prb(2)*Sig.^2+prb(3)*Sig+prb(4);
% meanLy=-aL.*gamma((bL+1)./bL); %%% the mean value of the Weibull distribution is a*Gamma((b+1)/b)
% meanRy=aR.*gamma((bR+1)./bR);
% 
% load('Zxya0d9epsilon1em4dt1em6.mat') %% Phase sensitivity function result obtained by the adjoint method. The adjoint method can be referred from e.g., Nakao 2016 Contemporary Physics and code therein.
% labelL=find(output(:,2)<0);
% theL=output(labelL,1);YtheL=output(labelL,3);
% meanL=zeros(1,length(Sig));
% for ii=1:length(Sig)
%     [~,lab]=min(abs(meanLy(ii)-YtheL));
%     meanL(ii)=theL(lab);
% end
% labelR=find(output(:,2)>0);
% theR=output(labelR,1);YtheR=output(labelR,3);
% meanR=zeros(1,length(Sig));
% for ii=1:length(Sig)
%     [~,lab]=min(abs(meanRy(ii)-YtheR));
%     meanR(ii)=theR(lab);
% end
% clear labelL labelR output theL theR YtheL YtheR
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% load('KtheSigma0d01and0d012.mat')
% load('Zxya0d9epsilon1em4dt1em6.mat')
% load('Ktheta(JumpStrength).mat')
% outputN=output(end:-1:1,:);
% %     Gammatheta=zeros(dN,1);
% Nz=zeros(1,2);
%     for ii=1:length(meanL)
%         thejumpL=meanL(ii);thejumpR=meanR(ii);
%         [~,labR]=min(abs(outputN(:,1)-thejumpR));
%         [~,labL]=min(abs(outputN(:,1)-thejumpL));
%         [~,labJR]=min(abs(thetar-thejumpR));[~,labRj]=min(abs(outputN(:,1)-thetaljump(labJR)));
%         [~,labJL]=min(abs(thetal-thejumpL));[~,labLj]=min(abs(outputN(:,1)-thetarjump(labJL)));
%         LabelSum=[1:1:labR labRj:1:labL labLj:1:length(outputN(:,1))];
%         magnifrac=2*pi/(2*pi-(outputN(labRj,1)-outputN(labR,1)+outputN(labLj,1)-outputN(labL,1)));
%         outputZ=[outputN(LabelSum,5)*magnifrac,outputN(LabelSum,1),outputN(LabelSum,3)];% Zy, phi, y
%         if ii==1
%             ZythetaN1=outputZ; 
%             Nz(ii)=size(outputZ,1);
%         else
%             ZythetaN2=outputZ;
%             Nz(ii)=size(outputZ,1);
%         end
%         
%     end
%     clear output outputN outputZ
%     
%     ZythetaN1(:,2)=linspace(0,2*pi,Nz(1));
%     ZythetaN2(:,2)=linspace(0,2*pi,Nz(2));
%     NGamma=1000;
%     thetaset=linspace(0,2*pi,NGamma);
%     ZythetaN=zeros(NGamma,4);
%     for ii=1:NGamma
%         [~,lab]=min(abs(ZythetaN1(:,2)'-thetaset(ii)));ZythetaN(ii,1:2)=ZythetaN1(lab,[1,3]); %%% Zy, y
%         [~,lab]=min(abs(ZythetaN2(:,2)'-thetaset(ii)));ZythetaN(ii,3:4)=ZythetaN2(lab,[1,3]);
%     end
%     ZythetaN=[ZythetaN;ZythetaN];%%% double period
%     Gammatheta=zeros(2,NGamma);
%     for ii=1:NGamma
%         Gammatheta(1,ii)=1/2/pi*sum(ZythetaN(ii:NGamma+ii-1,1).*(ZythetaN(1:NGamma,4)-ZythetaN(ii:NGamma+ii-1,2)))*2*pi/NGamma;
%         Gammatheta(2,ii)=1/2/pi*sum(ZythetaN(1:NGamma,3).*(ZythetaN(ii:NGamma+ii-1,2)-ZythetaN(1:NGamma,4)))*2*pi/NGamma;
%     end
%     Gammathetad=(Gammatheta(1,:)-Gammatheta(2,:));
%     plot(thetaset,Gammathetad)
%% (5) reduced and averaged model
% load('Gamma12andGammad.mat') %%% obtained from the above (4)
% load('KtheSigma0d01and0d012.mat') %% obtained from the code file 'ThetaJumpleft2right.m'
% omega1=omega/(1-Kthe1/2/pi);
% omega2=omega/(1-Kthe2/2/pi);
% delta=omega1-omega2;
% mu=[0.02,0.15];
% 
% muc=mu(2);
% phi0=6;dt=0.01;ts=0;te=200;
% [t,phi]=EulerPhi(phi0,dt,ts,te,Gammathetad,delta,muc);
% plot(t(1:10:end),phi(1:10:end),'-k')
% hold on

%% (6) original model
% ts=0;te=200;
% h=1e-5;
% load('Zxya0d9epsilon1em4dt1em6.mat')
% Sig=[0.01,0.012];
% mu=0.15;
% phi0=0;
% [~,lab2]=min(abs(output(:,1)-0));
% [~,lab1]=min(abs(output(:,1)-phi0));
% xy1=output(lab1,2:3);
% xy2=output(lab2,2:3);
% 
% 
% outputPositiveX=output(output(:,2)>=0,[1,3]);
% outputNegativeX=output(output(:,2)<0,[1,3]);
% outputPositiveX=outputPositiveX(1:50:end,:);
% outputNegativeX=outputNegativeX(1:50:end,:);
% 
% samplenum=100;xynDist=1000;
% phi=zeros(2*(floor((te-ts)/h/xynDist)+1),samplenum);
% parfor ss=1:samplenum
%     [~,phi(:,ss)]=par4originalphi(ss,xy1,xy2,h,ts,te,Sig,mu,xynDist,outputPositiveX,outputNegativeX);
% end
% t=linspace(ts,te,size(phi,1)/2);
% plot(t,phi(1:size(phi,1)/2,1))
% hold on
% plot(t,phi(size(phi,1)/2+1:end,1))

% % plot(t(1:1000:end),xyN(:,4))
% % hold on
% % plot(t(1:1000:end),xyN(:,2))


%% (7) reduced phase equation before averaging 
% load('LinearAndThirdOrderFitOfAandBversusSigma4Weibull.mat')
% Sig=[0.01,0.012];
% aL=pla(1)*Sig+pla(2);
% bL=plb(1)*Sig.^3+plb(2)*Sig.^2+plb(3)*Sig+plb(4);
% aR=pra(1)*Sig+pra(2);
% bR=prb(1)*Sig.^3+prb(2)*Sig.^2+prb(3)*Sig+prb(4);
% yy1=-0.667:0.0001:0;
% yy2=0:0.0001:0.667;
% CDFL=[1-exp(-(abs(yy1)/aL(1)).^bL(1));1-exp(-(abs(yy1)/aL(2)).^bL(2))];
% CDFR=[1-exp(-(yy2/aR(1)).^bR(1));1-exp(-(yy2/aR(2)).^bR(2))];
% load('Zxya0d9epsilon1em4dt1em6.mat')
% outputLTheta=output(output(:,2)<0,[1,3]);
% outputRTheta=output(output(:,2)>0&output(:,3)>0,[1,3]);
% Theleft=zeros(1,length(yy1));
% Theright=zeros(1,length(yy2));
% for ii=1:length(yy1)
%     [~,lab]=min(abs(yy1(ii)-outputLTheta(:,2)));
%     Theleft(ii)=outputLTheta(lab,1);
% end
% for ii=1:length(yy2)
%     [~,lab]=min(abs(yy2(ii)-outputRTheta(:,2)));
%     Theright(ii)=outputRTheta(lab,1);
% end
%%%%%%%%%%%%%% jump position CDF random generation
% load('CDFYandThetaSigma0d01and0d012dy0d0001.mat')
% load('Ktheta(JumpStrength).mat')
% Nsample=2e4;
% JumpLtheta=zeros(2,Nsample);
% JumpRtheta=zeros(2,Nsample);
% for kk=1:2
%     UsamL=rand(Nsample,1);
%     UsamR=rand(Nsample,1);
%     thetajumpL=zeros(Nsample,1);
%     thetajumpR=zeros(Nsample,1);
%     
%     for ii=1:Nsample
%         [~,lab]=min(abs(UsamL(ii)-CDFL(kk,:)));
%         thetajumpL(ii)=Theleft(lab);
%         [~,lab]=min(abs(UsamR(ii)-CDFR(kk,:)));
%         thetajumpR(ii)=Theright(lab);
%     end
%     JumpLtheta(kk,:)=thetajumpL;
%     JumpRtheta(kk,:)=thetajumpR;
% end
% load('Ktheta(JumpStrength).mat')
% thetakL=thetal(labL);thetakR=thetar(labR);
% JumpLtheta=JumpLtheta(:);JumpRtheta=JumpRtheta(:);
% JumpLthetaR=zeros(length(JumpLtheta),1);JumpRthetaL=zeros(length(JumpRtheta),1);
% for ii=1:length(JumpLtheta)
%     [~,lab]=min(abs(JumpLtheta(ii)-thetakL));
%     JumpLthetaR(ii)=thetakL(lab)+KthetaL(lab);
% end
% for ii=1:length(JumpRtheta)
%     [~,lab]=min(abs(JumpRtheta(ii)-thetakR));
%     JumpRthetaL(ii)=thetakR(lab)+KthetaR(lab);
% end
% JumpLtheta=reshape(JumpLtheta,2,[]);
% JumpRtheta=reshape(JumpRtheta,2,[]);
% JumpLthetaR=reshape(JumpLthetaR,2,[]);
% JumpRthetaL=reshape(JumpRthetaL,2,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% load('Zxya0d9epsilon1em4dt1em6.mat')
% outputTheandYandZy=output(end:-1:1,[1,2,3,5]);
% clearvars -except omega outputTheandYandZy thetakL thetakR KthetaL KthetaR
% load('JumpLRthetaSig0d01and0d012N20000.mat')
% 
% phi0=[0,0];
% h=0.01;ts=0;te=200;
% Sig=[0.01,0.012];
% mu=0.15;
% samplenum=100;
% t=ts:h:te;
% phi=zeros(2*length(t),samplenum);
% parfor jj=1:samplenum
%     [~,phi(:,jj)]=EulerCoupledReducedBeforeAverage(jj,samplenum,omega,phi0,h,ts,te,mu,outputTheandYandZy,JumpLtheta,JumpRtheta,JumpLthetaR,JumpRthetaL);
% end
% 
