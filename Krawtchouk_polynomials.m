function K=Krawtchouk_polynomials(n,N,p)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%   Downloaded from                                                            %
%   https://github.com/JSaulRivera/Computation-of-2D-and-3D-high-order-        %
%   discrete-orthogonal-moments.git                                            %
%                                                                              %
%                                                                              %
%   This code calculate the discrete orthogonal Krawtchouk polynomial          % 
%   for high order, using Gram-smith orthogonalization.'                       %
%                                                                              %
%                                                                              %
%                                                                              %
%   Please cite the following documents as:                                    %
%                                                                              %
%  *  José S. Rivera-Lopez, Cesar Camacho-Bello, and Lucia                     %
%     Gutiérrez-Lazcano, Chapter 3: "Computation of 2D and 3D High-order       %
%     Discrete Orthogonal Moments". Recent Progress in Image Moments and       %
%     Moment Invariants, GCSR Volume 7 (2021), 53-74, DOI: 10.15579/gcsr.      %
%     vol7.ch3.                                                                %
%                                                                              %
%        Bibtex:                                                               %
%        @article{Rivera2021Recent,                                            %
%        title={Recent Progress in Image Moments and Moment Invariants},       %
%        author={José S. Rivera-Lopez and César Camacho-Bello and Lucia        %
%        Gutiérrez-Lazcano},                                                   %
%        journal={Science Gate},                                               %
%        volume={7},                                                           %
%        pages={53--74},                                                       %
%        year={2021},                                                          %
%        URL={https://sciencegatepub.com/sgp-books/gcsr/gcsr_vol7/},           %
%        ISBN={2241-9063}                                                      %
%        }                                                                     %
%                                                                              %
%                                                                              %
%  *  C. Camacho-Bello and J. S. Rivera-Lopez, â€œSome computational aspects   %
%     of tchebichef moments for higher orders,â€?Pattern Recognition           %
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

x=1:N;
n=n-1;
w = binopdf(x,N,p);
   

K(1,:)= sqrt(w);
K(2,:)= ((x-p*N).* K(1,:)/sqrt(N*p*(1-p)));
   

i=1;
w1=sqrt((N-i+1)*i);


for i=2:n
    
    w2=sqrt((N-i+1)*i);
    wn=(x-i+1-p*(N-2*i+2))/sqrt(p*(1-p));
   
    K(i+1,:)=(wn.*K(i,:)-w1.*K(i-1,:))/w2;
    An=An2     ;
      
end

h=sqrt(sum(K(n+1,:).^2))+eps;
K2(n+1,:)=K(n+1,:);
  
     
for k=1:n
    K(n+1,:)=K(n+1,:)-(K2(n+1,:)*K(k,:)').*K(k,:);  
end


K(n+1,:)=K(n+1,:)./h;
end
