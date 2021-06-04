function C=Charlier_polynomials(n,N,a)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%   Downloaded from                                                            %
%   https://github.com/JSaulRivera/Computation-of-2D-and-3D-high-order-        %
%   discrete-orthogonal-moments.git                                            %
%                                                                              %
%                                                                              %
%   This code calculate the discrete orthogonal Charlier polynomial            % 
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
   n=n-1;
   x=0:N-1;
   
   C(1,1)=sqrt(exp(-a));
   
   for i=2:N
  
       C(1,i)=C(1,i-1)*sqrt(a/(i-1)); 
   end
   
   C(2,:)= ((a-x)/sqrt(a)).*C(1,:);


   h=sqrt(sum(C(1,:).^2))+eps;
   C(1,:)=C(1,:)./h;
   
   h=sqrt(sum(C(2,:).^2))+eps;
   C2(2,:)=C(2,:);
   
   C(2,:)=C(2,:)-(C2(2,:)*C(1,:)').*C(1,:);  
   C(2,:)=C(2,:)./h;
   
   An=1;
     
   for i=2:n
    
      An1=sqrt(i);
      Bn=(a-(x-1)+i-1)/sqrt(a);
   
      C(i+1,:)=(Bn.*C(i,:)-An.*C(i-1,:))/An1;
      An=An1;  
      
      h=sqrt(sum(C(i+1,:).^2))+eps;
      C2(i+1,:)=C(i+1,:);
   
      
         for k=1:i
             C(i+1,:)=C(i+1,:)-(C2(i+1,:)*C(k,:)').*C(k,:);  
         end
      
      C(i+1,:)=C(i+1,:)./h;
   end
   
   

end
