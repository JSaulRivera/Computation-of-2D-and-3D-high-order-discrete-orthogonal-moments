function H=Hahn_polynomials(n,M,a,b)
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%   Downloaded from                                                            %
%   https://github.com/JSaulRivera/Computation-of-2D-and-3D-high-order-        %
%   discrete-orthogonal-moments.git                                            %
%                                                                              %
%                                                                              %
%   This code calculate the discrete orthogonal Hahn polynomial                % 
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

    n=n-1;
    for x=1:M
        H(1,x)= sqrt(HahnW(a,b,M,x-1))/(Hahndn2(a,b,M,0));
        H(2,x)= ((a+b+2)*(x-1)-(b+1)*(M-1))*sqrt(HahnW(a,b,M,x-1))/...
            (Hahndn2(a,b,M,1));  
    end


    q=n;


    An= sqrt((a+b+1)*(M-1)*(a+1)*(b+1)*(a+b+1+M))/((a+b+2)*sqrt(a+b+1));
    
    
    for n=2:q

        An1=(sqrt((n-1)*(a+n-1)*(b+n-1)*(a+b+n-1)*(M-n+1)*(a+b+(2*n)+1)))...
            /((a+b+(2*n)-2));

        for x=0:M-1
      
        Bn(x+1)=sqrt(a+b+2*n+1)*(x-((a-b+2*M-2)/4)-(((b^2-a^2)*(a+b+2*M))/(4*(a+b+2*n-2))));
          
        end

        H(n+1,:)=(Bn.*H(n,:)-An1*H(n-1,:))/An; 
 
        An=An1;
        
        h=sqrt(sum(H(n+1,:).^2))+eps;
        H2(n+1,:)=H(n+1,:);
  
   
   
       for k=1:n
            H(n+1,:)=H(n+1,:)-(H2(n+1,:)*H(k,:)').*H(k,:);  
       end
 
   
       H(n+1,:)=H(n+1,:)./h;
   
       end


end

function W=HahnW(alpha,beta,N,x)
   W=Pochhammer(x+1,beta)*Pochhammer(N-x,alpha);
end

function P= Pochhammer(a,k)
c=a;
if k>1
   for h=1:k-1
     c=c*(a+h);
     %disp(c);
   end
end

P=c;
end


function D=Hahndn2(a,b,N,n)
%D=(Pochhammer((n+1),b)/Pochhammer((a+n+1),b))*(Pochhammer(N-n,a+b+2*n+1)/(a+b+2*n+1));
% D=(Pochhammer((n+1),b)*(Pochhammer2(N-n,a+b+2*n+1,a+n+1,b))/(a+b+2*n+1));
D=sqrt(Pochhammer((n+1),b)/(a+b+2*n+1))*Pochhammer2(N-n,a+b+2*n+1,a+n+1,b);
end

function P= Pochhammer2(a,k,a2,k2)
% c=a;
% c2=a2;
c3=sqrt(a/a2);
if k>1
   for h=1:k-1
       
     if h<=k2-1
       c3=c3*sqrt((a+h)/(a2+h));
       
     else
       c3=c3*sqrt(a+h);
     end
     
     %disp(c);
   end
end

P=c3;
end