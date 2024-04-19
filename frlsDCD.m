function [w,erro] = frlsDCD(lambda,H,u,d,M,delta,Nu)
%#codegen
N = length(u);
w = zeros(M,1);
beta = zeros(M,1);
x = zeros(M,1);
erro = zeros(1,N);

R=[delta*ones(1,M);zeros(M-1,M)];

indri=ones(M,1);

indrmod=M:-1:1;

for n=1:N
    x=[u(n);x(1:M-1)];
    erro(n)=d(n) - w'*x;
    for i=1:M
        temp=indri(i);
        indri(i)= mod1(indri(i)-1,indrmod(i)); 
        R(i,indri(i))=lambda*R(i,temp) + x(1)*x(i);
    end
   beta = lambda*beta + erro(n)*x;
   [dw,beta]=fdcd(M,H,16,Nu,R,indri,indrmod,beta);
   w=w+dw;
end
end

function [dw,res]=fdcd(M,H,B,Nu,R,indri,indrmod,beta)
    h=H/2;
    b=1;
    res=beta;
    dw=zeros(M,1);
    for k=1:Nu
        [~,p] = max(abs(res));
        pp=mod1(indri(1)+p-1,indrmod(1));
        while abs(res(p)) <= (h/2)*R(1,pp)
            b=b+1;
            h=h/2;
            if b> B
                break
            end
        end
        dw(p)=dw(p)+sign(res(p))*h;
        fact=sign(res(p))*h;
        j=p;
        dir=-1;
        k=0;
        kstep=1;
        for i=1:M
            res(i)=res(i)-fact*R(j,mod1(indri(j)+k,indrmod(j)));
            j=j+dir;
            k=k+kstep;
            if j==0
                j=2;
                dir=+1;
                k=p-1;
                kstep=0;
            end             
        end
    end
end

function y=mod1(x,m)
    
    if (x==m || x==0)
        y=m;
    else
        y=mod(x,m);
    end
end