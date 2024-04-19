function [w,erro]=nlms(mu,u,d,M,delta)
%#codegen
% Normalized LMS
% Call:
% [W,erro,ea]=nlms(mu,u,d,M,delta);
%
% Argumentos de entrada:
% mu = passo de adaptação
% M = tamanho do filtro
% u = sinal de entrada
% delta = constante
%
% Argumentos de saida:
% erro = erro de estimação
% w = coeficientes finais do filtro
u=u(:);
d=d(:);
N = length(u);
U = zeros(M,1);
erro=zeros(1,N);
w=zeros(M,1);

for n=1:N 
    U=[u(n);U(1:M-1)];
    erro(n)=d(n)-w'*U;
    w = w + (mu/(delta+U'*U))*erro(n)*U;
end