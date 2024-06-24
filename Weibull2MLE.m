
%% %%%%%%%%%%% Left position %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramsuml=zeros(5,2);
for ii=1:5
% ii=1;
load(strcat('ypositionsigma0d0',num2str(ii),'.mat')) %%% obtained from Monte Carlo simulations
yleftsum=yleft(:);
yleftsum(yleftsum>-0.3)=[];
data=-yleftsum(:);
weibull2pdf = @(x,a,b) (b/a).*((x/a).^(b-1)).*exp(-(x/a).^b);%%% two parameter weibull distribution
opt = statset('MaxIter',1e5,'MaxFunEvals',1e5,'FunValCheck','off');
% 'MaxIter',1e5 — Increase the maximum number of iterations to 1e5.
%'MaxFunEvals',1e5 — Increase the maximum number of object function evaluations to 1e5.
% 'FunValCheck','off' — Turn off checking for invalid object function values.

params = mle(data,'pdf',weibull2pdf,'start',[0.5 10],...
     'Options',opt,'LowerBound',[0 0],'UpperBound',[Inf Inf]);

a1=params(1);b1=params(2);
dx=(max(data)-min(data))/200;
xx=min(data)-0.01:dx:max(data)+0.01;
estimateWeibull=b1/a1*(xx/a1).^(b1-1).*exp(-(xx/a1).^b1);
plot(-xx,estimateWeibull)
hold on
h1=histogram(-data);
h1.Normalization='pdf';
paramsuml(ii,:)=params;
end
%% %%%%%%%%%%% Right position %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramsumr=zeros(5,2);
for ii=1:5
% ii=1;
load(strcat('ypositionsigma0d0',num2str(ii),'.mat')) %%% obtained from Monte Carlo simulations
yrightsum=yright(:);
yrightsum(yrightsum<0.3)=[];
data=yrightsum(:);
weibull2pdf = @(x,a,b) (b/a).*((x/a).^(b-1)).*exp(-(x/a).^b);%%% two parameter weibull distribution
opt = statset('MaxIter',1e5,'MaxFunEvals',1e5,'FunValCheck','on');
% 'MaxIter',1e5 — Increase the maximum number of iterations to 1e5.
%'MaxFunEvals',1e5 — Increase the maximum number of object function evaluations to 1e5.
% 'FunValCheck','off' — Turn off checking for invalid object function values.

params = mle(data,'pdf',weibull2pdf,'start',[0.1 1.5],...
     'Options',opt,'LowerBound',[0 0],'UpperBound',[Inf Inf]);

a1=params(1);b1=params(2);
dx=(max(data)-min(data))/200;
xx=min(data)-0.01:dx:max(data)+0.01;
estimateWeibull=b1/a1*(xx/a1).^(b1-1).*exp(-(xx/a1).^b1);
plot(xx,estimateWeibull)
hold on
h1=histogram(data);
h1.Normalization='pdf';
paramsumr(ii,:)=params;
end
