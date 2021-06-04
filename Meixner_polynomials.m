function M=Meixner_polynomials(n,N,u,b)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%   Downloaded from                                                            %
%   https://github.com/JSaulRivera/Computation-of-2D-and-3D-high-order-        %
%   discrete-orthogonal-moments.git                                            %
%                                                                              %
%                                                                              %
%   This code calculate the discrete orthogonal Meixner polynomial             % 
%   for high order, using Gram-smith orthogonalization.'                       %
%                                                                              %
%                                                                              %
%                                                                              %
%   Please cite the following documents as:                                    %
%                                                                              %
%  *  José S. Rivera-Lopez, César Camacho-Bello, and Lucia                     %
%     Gutiérrez-Lazcano, Chapter 3: "Computation of 2D and 3D High-order       %
%     Discrete Orthogonal Moments". Recent Progress in Image Moments and       %
%     Moment Invariants, GCSR Volume 7 (2021), 53-74, DOI: 10.15579/gcsr.      %
%     vol7.ch3.                                                                %
%                                                                              %
%        Bibtex:                                                               %
%        @article{Rivera2021Recent,                                            %
%        title={Recent Progress in Image Moments and Moment Invariants},       %
%        author={JosÃ© S. Rivera-Lopez and CÃ©sar Camacho-Bello and Lucia      %
%        GutiÃ©rrez-Lazcano},                                                  %
%        journal={Science Gate},                                               %
%        volume={7},                                                           %
%        pages={53--74},                                                       %
%        year={2021},                                                          %
%        URL={https://sciencegatepub.com/sgp-books/gcsr/gcsr_vol7/},           %
%        ISBN={2241-9063}                                                      %
%        }                                                                     %
%                                                                              %
%                                                                              %
%  *  C. Camacho-Bello and J. S. Rivera-Lopez, Some computational aspects      %
%     of tchebichef moments for higher orders, Pattern Recognition             %
%     Letters, vol. 112, pp. 332â€“339, 2018.                                  %
%                                                                              %
%        Bibtex:                                                               % 
%        @article{camacho2018some,                                             %
%        title={Some computational aspects of Tchebichef moments for           %
%        higher orders},                                                       %
%        author={Camacho-Bello, C{\'e}sar and Rivera-Lopez, Jos{\'e} S},       %
%        journal={Pattern Recognition Letters},                                %
%        volume={112},                                                         %
%        pages={332--339},                                                     %
%        year={2018},                                                          %
%        publisher={Elsevier}                                                  %
%        }                                                                     %                                                    
%                                                                              %
%                                                                              %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=n-1;


for x=1:N   
 M(1,x)= sqrt(Pochhamer(x,b-1)*(((u^(x-1)))/(factorial(b-1)))*((1-u)^b));
 M(2,x)= (b+(x-1)-((x-1)/u))*sqrt(M(1,x)/(b));
end


h=sqrt(sum(M(1,:).^2))+eps;
M(1,:)=M(1,:)./h;
   
h=sqrt(sum(M(2,:).^2))+eps;
M2(2,:)=M(2,:);
   
M(2,:)=M(2,:)-(M2(2,:)*M(1,:)').*M(1,:);  
M(2,:)=M(2,:)./h;
   
 
q=n;
w1=sqrt(b);  

for n=2:q
   
   w2=sqrt((n)*(b+n-1));
    
   for x=1:N 
        
       w =- (((x-1)-((x-1)*u)-n+1-(u*n)+u-(b*u))*sqrt(u))/u;
       M(n+1,x)=(w.*M(n,x)-w1*M(n-1,x))/w2;

   end
    
   w1=w2;  
    
   h=sqrt(sum(M(n+1,:).^2))+eps;
   M2(n+1,:)=M(n+1,:);
  
   
   
   for k=1:n
      M(n+1,:)=M(n+1,:)-(M2(n+1,:)*M(k,:)').*M(k,:);  
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
