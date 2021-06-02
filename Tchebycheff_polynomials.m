function T=PolyTchebichef(q,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%   Downloaded from                                                            %
%   https://github.com/JSaulRivera/Computation-of-2D-and-3D-high-order-        %
%   discrete-orthogonal-moments.git                                            %
%                                                                              %
%                                                                              %
%   This code calculate the discrete orthogonal Tchebichef polynomial          % 
%   for high order, using Gram-smith orthogonalization.'                       %
%                                                                              %
%                                                                              %
%                                                                              %
%   Please cite the following documents as:                                    %
%                                                                              %
%  *  José S. Rivera-Lopez, César Camacho-Bello, and Lucia                     %
%     Gutiérrez-Lazcano, Chapter 3: “Computation of 2D and 3D High-order       %
%     Discrete Orthogonal Moments”. Recent Progress in Image Moments and       %
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
%  *  C. Camacho-Bello and J. S. Rivera-Lopez, “Some computational aspects     %
%     of tchebichef moments for higher orders,”Pattern Recognition             %
%     Letters, vol. 112, pp. 332–339, 2018.                                    %
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

q=q-1;
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
