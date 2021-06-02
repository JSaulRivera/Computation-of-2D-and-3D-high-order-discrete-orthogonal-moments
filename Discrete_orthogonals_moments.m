

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%   Downloaded from                                                            %
%   https://github.com/JSaulRivera/Computation-of-2D-and-3D-high-order-        %
%   discrete-orthogonal-moments.git                                            %
%                                                                              %
%                                                                              %
%   This code calculate the discrete orthogonal Moments for high order         % 
%                                                                              %
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

clc
clearvars
close all

%Read image
imagen= imread('Test_image.jpeg');
[N M]=size(imagen(:,:,1));

% Order moments
 n=N;

 % Obtain Polynomial base
 Polynomial_base = Tchebycheff_polynomials(n,N);
    
 Moments=zeros(n,n,3);
    
    for i=1:3
        I=double(imagen(:,:,i));
        % computation of discrete orthogonal moments
        Moments(:,:,i)=Polynomial_base*I*Polynomial_base';
        %  calculation for image reconstruction with discrete orthogonal moments
        Image_reconstructed(:,:,i)=Polynomial_base'*Moments(:,:,i)*Polynomial_base;
    end
    
 % Mostrar Imagen reconstruida
 figure(2)
 imshow(uint8(Image_reconstructed))

clc
clearvars
