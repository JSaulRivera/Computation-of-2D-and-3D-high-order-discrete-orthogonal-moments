function M=Meixner_polynomials(n,N,a1,b)
n=n-1;


for x=1:N
    
 M(1,x)= sqrt(Pochhamer(x,b-1)*(((a1^(x-1)))/(factorial(b-1)))*((1-a1)^b));
 M(2,x)= (b+(x-1)-((x-1)/a1))*sqrt(M(1,x)/(b));
end


   h=sqrt(sum(M(1,:).^2))+eps;
   M(1,:)=M(1,:)./h;
   
   h=sqrt(sum(M(2,:).^2))+eps;
   M2(2,:)=M(2,:);
   
   M(2,:)=M(2,:)-(M2(2,:)*M(1,:)').*M(1,:);  
   M(2,:)=M(2,:)./h;
   
 
 q=n;
  An=sqrt(b);  
for n=2:q
    

    
    An2=sqrt((n)*(b+n-1));
      for x=1:N 
%        
   BD =- (((x-1)-((x-1)*a1)-n+1-(a1*n)+a1-(b*a1))*sqrt(a1))/a1;
   M(n+1,x)=(BD.*M(n,x)-An*M(n-1,x))/An2;

      end
    An=An2;  
    
  
    
    
    
    
    
   h=sqrt(sum(M(n+1,:).^2))+eps;
   M2(n+1,:)=M(n+1,:);
  
   
   if n>1
       for k=1:n
         M(n+1,:)=M(n+1,:)-(M2(n+1,:)*M(k,:)').*M(k,:);  
       end
   end
   
   M(n+1,:)=M(n+1,:)./h;
   
end

end

function [p] = Pochhamer(num, N)
p=1;
for i=0:N-1
    p=p*(num+i);
end
end