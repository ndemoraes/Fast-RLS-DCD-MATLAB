function [w,erro] = rlsDCD(lambda,H,u,d,M,delta,Nu)
%#codegen
N = length(u);
w = zeros(M,1);
beta = zeros(M,1);
x = zeros(M,1);
erro = zeros(1,N);
R = delta*eye(M);
for n=1:N
    x = [u(n); x(1:M-1)];
    erro(n) = d(n) - w'*x;
    % dessa maneira só é feita a atualização da primeira coluna da matriz R,
    % aproveitando a sua simetria.
    for j=M-1:-1:1
        for i=1:M-1
            R(i+1,j+1) = R(i,j);
        end
        R(j+1,1) = lambda*R(j+1,1) + x(1) * x(j+1);
        R(1,j+1) = R(j+1,1);
    end
    R(1,1) = lambda * R(1,1) + x(1)^2;
    %%%%%%%%%
    beta = lambda*beta + erro(n)*x;
    [dw,beta] = DCD(M,beta,R,H,Nu);
    w = w + dw;
end
end

function [dw,beta] = DCD(M,beta,R,H,Nu)
%%%%%% DCD %%%%%%
dw_res = zeros(M,1);
res=beta;
h=H/2;
b=1;
B=16;
for k=1:Nu
    [a,p] = max(abs(res));
    while a <= (h/2)*R(p,p)
        b=b+1;
        h=h/2;
        if b>B
            break
        end
    end
    dr = sign(res(p))*h;
    dw_res(p) = dw_res(p) + dr;
    res = res - dr*R(:,p);
end
beta=res;
dw=dw_res;
end