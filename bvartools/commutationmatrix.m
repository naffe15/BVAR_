function K = commutationmatrix(n)

K = zeros(n^2,n^2);

for n1 = 1 : n
    for n2 = 1 : n
       tmp_ = zeros(n);
       tmp_(n2,n1) = 1;
       y1 = (n1-1)*n+1 : n1*n;
       y2 = (n2-1)*n+1 : n2*n;
       K(y1, y2) = tmp_;       
    end
end



end