function D = duplicationmatrix(n)
% Duplication matrix: vec(P)=Dvech(P)
% tic
m=1/2*n*(n+1);
nsq=n^2;
DT=sparse(m,nsq);
for j=1:n
    for i=j:n
        ijth=(j-1)*n+i;
        jith=(i-1)*n+j;
        vecTij=sparse(ijth,1,1,nsq,1);
        vecTij(jith,1)=1;
        k=(j-1)*n+i-1/2*j*(j-1);
        uij=sparse(k,1,1,m,1);
        DT=DT+uij*vecTij';
    end
end
% D=DT';
D = full(DT');