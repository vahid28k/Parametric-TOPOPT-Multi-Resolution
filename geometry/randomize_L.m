function [iter_U]=randomize_L(nelx,nely,penal,rmin,inf,x)

%% INPUT
% nelx=100;
% nely=40;
% volfrac=0.5;
% penal=3;
% rmin=2;
% inf=0;
% theta=0;
mag=1;
Ey=1;
ft=1;

%% MATERIAL PROPERTIES
E0 = 1*Ey;

Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
dx=1;   %size of the element in x
dy=1;   %size of the element in y
%stiffness matrix of one element (rectangular)
KE=stiff_ele(E0,nu,dx,dy);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

fixeddofs=[1 2];
for i=1:nelx
    fixeddofs=[fixeddofs 1+i*(nely+1)*2 2+i*(nely+1)*2];
end;
%fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION

xPhys = x;
loop = 0;
change = 1;
loop = loop + 1;
  %% FE-ANALYSIS
  xfilt=(H*reshape(xPhys,nelx*nely,1))./Hs;
  
 load data_100.mat;load des_quad.mat;wquad=w';  
 if nelx==4 out=out_4; end;
 if nelx==10 out=out_10; end;
 if nelx==20 out=out_20; end;
 if nelx==50 out=out_50; end;
 
  m=size(wquad,2);
  beta=2;eta=out*0.1+0.45;
  for i=1:m
      
      for j=1:(nelx*nely)
         xthr(j,1) = (tanh(beta*eta(j,i)) + tanh(beta*(xfilt(j) - eta(j,i))) ) /  (tanh(beta*eta(j,i)) + tanh(beta*(1 - eta(j,i))) );
         deri_x(j,1) = ( beta * (sech(beta*(xfilt(j)-eta(j,i))))^2 ) / (tanh(beta*eta(j,i)) + tanh(beta*(1 - eta(j,i))) ); 
       end;
       xPhys=reshape(xthr(:,1),nely,nelx);
      
      F = sparse(2*(nely+1)*(nelx+1),1);
      F(2*(nelx+1)*(nely+1)-nely+inf,1) = 1;%ver
      F(2*(nelx+1)*(nely+1)-nely+inf-1,1) = 0.2;%hor
      U = zeros(2*(nely+1)*(nelx+1),1);
      
      sK = reshape(KE(:)*( [Emin+(xPhys(:)'.^penal).*(E0-Emin)]),64*nelx*nely,1);%%%%%%%%%%%%%
      K = sparse(iK,jK,sK); K = (K+K')/2;
      U(freedofs) = K(freedofs,freedofs)\F(freedofs);
      
      %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS

       
      iter_U(:,i)=U;

      
  end;
  
end
 
%element matrix function
function KE=stiff_ele(E,nu,dx,dy)
%Gauss Points (4 in total)
egv=[-1/sqrt(3),1/sqrt(3)]; %2 gauss point x
ngv=egv;                    %2 gauss point y
wg=[1,1];                   %weight

%Constitutive matrix
matC=consti(E,nu);
%initialize stiffness matrix
KE=zeros(8,8);
%Gauss quadrature integration
for eit=1:2
    eg=egv(eit);
    for nit=1:2
        ng=ngv(nit);
        %Derivative of the shape functions
        [Be,Jdet]=Der_shape_fun(eg,ng,dx,dy);
        %Derivative of the shape function for displacement in x and y
        BKe=[Be(1,1) 0 Be(2,1) 0 Be(3,1) 0 Be(4,1) 0;...
               0 Be(1,1) 0 Be(2,1) 0 Be(3,1) 0 Be(4,1);
                Be(1,2) 0 Be(2,2) 0 Be(3,2) 0 Be(4,2) 0;...
                0 Be(1,2) 0 Be(2,2) 0 Be(3,2) 0 Be(4,2)];
        %Stifness for element in domain 1
        KE=KE+wg(eit)*wg(nit)*BKe'*matC*BKe*Jdet;
    end
end
end

%Derivative of Shape function for [dN1/dx dN2/dx ...;dN1/dy dN2/dy ... ]
%and determinant of Jacobian
function [Be,Jdet]=Der_shape_fun(eg,ng,dx,dy)
%coordinates of the 4 nodes (specific case for rectangle)
Coorde=[0 dx dx 0;0 0 dy dy]; 

%Derivative of the shape function
DNe=1/4*[-(1-ng) 1-ng 1+ng -(1+ng);...
    -(1-eg) -(1+eg) 1+eg 1-eg]';

%Jacobian of the element e at location e,n
J=Coorde*DNe;
Jdet=det(J);

%Derivative of the shape function
Be=DNe/J;
end


%Constitutive relationship in matrix form
function matC=consti(E,nu)
%Lame constants
lam=E*nu/(1+nu)/(1-2*nu);   
miu=E/(2*(1+nu));
lam=2*miu*lam/(lam+2*miu);  %plane stress
%Constitutive matrix
matC=[lam+2*miu 0 0 lam;0 miu miu 0;0 miu miu 0;lam 0 0 lam+2*miu];
end

