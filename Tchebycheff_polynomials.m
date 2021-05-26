function T=PolyTchebichef(q,N)

x=0:N-1;
w=2*x-N+1;
n=1;
w1=n*sqrt((N^2-n^2)/((2*n+1)*(2*n-1)));

T(1,:)=ones(size(x,2),1).*(1/sqrt(N));
T(2,:)=(w./w1).*T(1,:);


for n=2:q
   
    w2=n*sqrt((N^2-n^2)/((2*n+1)*(2*n-1)));
    T(n+1,:)=(w./w2).*T(n,:)-(w1./w2).*T(n-1,:);
    w1=w2;

   T2=T(n+1,:);
  
       for k=1:n
         T(n+1,:)=T(n+1,:)-(T2*T(k,:)').*T(k,:);  
       end

   h=sqrt(sum(T(n+1,:).^2))+eps;
   T(n+1,:)=T(n+1,:)./h;
   
end


end