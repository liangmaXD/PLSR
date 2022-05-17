function [B,Wstar,T,P,Q,W,R2X,R2Y]=plsnipals(X,Y,A)
%+++ The NIPALS algorithm for both PLS-1 (a single y) and PLS-2 (multiple Y)
%+++ X: n x p matrix
%+++ Y: n x m matrix
%+++ A: number of latent variables
%+++ Code: Hongdong Li, lhdcsu@gmail.com, Feb, 2014
%+++ reference: Wold, S., M. Sj鰏tr鰉, and L. Eriksson, 2001. PLS-regression: a basic tool of chemometrics,
%               Chemometr. Intell. Lab. 58(2001)109-130.



varX=sum(sum(X.^2));% X的平方累加
varY=sum(sum(Y.^2));% Y的平方累加
for i=1:A
    error=1;
    u=Y(:,1);
    niter=0;
    while (error>1e-8 && niter<1000)  % for convergence test 收敛性检验
        w=X'*u/(u'*u); 
        w=w/norm(w); %计算特征向量wi
        t=X*w; %计算成分ti的得分
        q=Y'*t/(t'*t);  % regress Y against t; %计算βi
        u1=Y*q/(q'*q);  % 
        error=norm(u1-u)/norm(u);
        u=u1;
        niter=niter+1;
    end
    p=X'*t/(t'*t);%计算αi
    X=X-t*p'; % 残差阵 X-X^=X-ti*αi
    Y=Y-t*q'; % 残差阵 Y-Y^=Y-ti*βi
    
    %+++ store
    W(:,i)=w;
    T(:,i)=t;
    P(:,i)=p;
    Q(:,i)=q;
    
end

%+++ calculate explained variance
R2X=diag(T'*T*P'*P)/varX;
R2Y=diag(T'*T*Q'*Q)/varY;

Wstar=W*(P'*W)^(-1); 
B=Wstar*Q';
Q=Q';

%+++ 
