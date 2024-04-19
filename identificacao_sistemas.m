clear; close all; clc;

epsilon=1.0;
delta_v=0.1;
lambda = 0.99;
delta = 1.0;

L=100;

M = [10, 100, 1000];
N=5000;
NM=length(M);
TRLS=zeros(NM,L);
TRLSDCD=zeros(NM,L);
TRLSDCD4=zeros(NM,L);
TfRLSDCD=zeros(NM,L);
TfRLSDCD4=zeros(NM,L);

for k=1:length(M)
    hi=randn(M(k),1);
    u=randn(N,1);
    d=filter(hi,1,u)+delta_v*randn(N,1);
    disp(['M= ', num2str(M(k))]);
    for i=1:L
        t2=@() rls(lambda,u,d,M(k),delta);
        t3=@() rlsDCD(lambda,4,u,d,M(k),delta,1);
        t4=@() rlsDCD(lambda,4,u,d,M(k),delta,4);
        t6=@() frlsDCD(lambda,4,u,d,M(k),delta,1);
        t7=@() frlsDCD(lambda,4,u,d,M(k),delta,4);
        
        TRLS(k,i)=timeit(t2);
        TRLSDCD(k,i)=timeit(t3);
        TRLSDCD4(k,i)=timeit(t4);
        TfRLSDCD(k,i)=timeit(t6);
        TfRLSDCD4(k,i)=timeit(t7);
        
        disp(['Tempo RLS                ', num2str(TRLS(k,i))]);
        disp(['Tempo RLSDCD Nu =1       ', num2str(TRLSDCD(k,i))]);
        disp(['Tempo RLSDCD Nu =4       ', num2str(TRLSDCD4(k,i))]);
        disp(['Tempo fRLSDCD Nu =1      ', num2str(TfRLSDCD(k,i))]);
        disp(['Tempo fRLSDCD Nu =4      ', num2str(TfRLSDCD4(k,i))]);
        
%         disp(['Tempo NLMS               ', num2str(MeNLMS(k,i))]);
%         disp(['Tempo RLS                ', num2str(MeRLS(k,i))]);
%         disp(['Tempo RLSDCD Nu =1       ', num2str(MeRLSDCD(k,i))]);
%         disp(['Tempo RLSDCD Nu =4       ', num2str(MeRLSDCD4(k,i))]);
%         disp(['Tempo RLSDCD Nu =16      ', num2str(MeRLSDCD16(k,i))]);
%         disp(['Tempo fRLSDCD Nu =1      ', num2str(MefRLSDCD(k,i))]);
%         disp(['Tempo fRLSDCD Nu =4      ', num2str(MefRLSDCD4(k,i))]);
%         disp(['Tempo fRLSDCD Nu =16     ', num2str(MefRLSDCD16(k,i))]);
    end
    disp('  ')
end
MeRLS=mean(TRLS,2);
MeRLSDCD=mean(TRLSDCD,2);
MeRLSDCD4=mean(TRLSDCD4,2);
MefRLSDCD=mean(TfRLSDCD,2);
MefRLSDCD4=mean(TfRLSDCD4,2);
MiRLS=mean(TRLS,2);
MiRLSDCD=min(TRLSDCD,2);
MiRLSDCD4=min(TRLSDCD4,2);
MifRLSDCD=min(TfRLSDCD,2);
MifRLSDCD4=min(TfRLSDCD4,2);

% save('fRLSDCD_mex','M','L','N','MefRLSDCD','MeRLSDCD','MefRLSDCD4','MeRLSDCD4','MefRLSDCD16','MeRLSDCD16','MeNLMS','MeRLS');

clf()
loglog(M,MeRLS)
hold on
loglog(M,MeRLSDCD)
loglog(M,MeRLSDCD4)
loglog(M,MefRLSDCD)
loglog(M,MefRLSDCD4)
hold off
legend('RLS','RLS-DCD Nu = 1','RLS-DCD Nu = 4',...
        'fast RLS-DCD Nu = 1','fast RLS-DCD Nu = 4')
grid()
xlabel("M")
title("Running time")



