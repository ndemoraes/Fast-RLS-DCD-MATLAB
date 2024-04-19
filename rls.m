function [w,erro] = rls(lambda,u,d,M,delta)
%#codegen
% Algoritmo RLS com simetria na atualizacao da matriz R
% Call:
% [w,erro] = rls(lambda,u,d,M,delta);
% 
% Argumentos de entrada:
% lambda = fator de esquecimento
% u = sinal de entrada do filtro
% d = sinal desejado
% M = tamanho do filtro.
% delta = fator de regularizacao
% 
% Argumentos de saida:
% erro = erro a priori
% w = caminho de adaptacao dos coeficientes do filtro
u=u(:);
d=d(:);
N = length(u);
w = zeros(M,1);
x = zeros(1,M);
P = (1/delta)*eye(M,M);
erro = zeros(1,N);
lb_inv = (1/lambda);
for n=1:N
    x = [u(n) x(1:M-1)];
	erro(n) = d(n) - x*w; % erro a priori
    g = P*x';
    gamma = 1/(lambda + x*g);
    k = gamma*g;
    w = w + k*erro(n);
    P = (lb_inv)*(P - (g*g')*gamma); % o parentesis em (g*g') garante a simetria, ou fazer g*g'*gamma.
end
end