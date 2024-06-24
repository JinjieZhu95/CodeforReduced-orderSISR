clear

%% correspondence of jumping position for theta
load('Zxya0d9epsilon1em4dt1em6.mat')%% Phase sensitivity function result obtained by the adjoint method. The adjoint method can be referred from e.g., Nakao 2016 Contemporary Physics and code therein.

labelxl=find(output(:,2)<0);
labelxr=find(output(:,2)>0);
yl=output(labelxl,3);yr=output(labelxr,3);
thetal=output(labelxl,1);thetar=output(labelxr,1);

% thetarjump=zeros(length(yl),1);
% for ii=1:length(yl)
%     [~,lab]=min(abs(yl(ii)-yr));
%     thetarjump(ii)=thetar(lab);
% end

thetaljump=zeros(length(yr),1);
for ii=1:length(yr)
    [~,lab]=min(abs(yr(ii)-yl));
    thetaljump(ii)=thetal(lab);
end

%%
labL=find(yl<0);
KthetaL=thetarjump(labL)-thetal(labL);

labR=find(yr>0);
KthetaR=thetaljump(labR)-thetar(labR);
plot(thetal(labL),KthetaL)
hold on
plot(thetar(labR),KthetaR)