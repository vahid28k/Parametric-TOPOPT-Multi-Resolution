function [iter_cdc]=randomize_L(nelx,nely,penal,x)

%% INPUT

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

fixeddofs=[2*nely+1 2*nely+2];
for i=1:nelx
    fixeddofs=[fixeddofs (2*nely+1)+[(nely+1)*2]*i (2*nely+2)+[(nely+1)*2]*i];
end;
    
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

%% INITIALIZE ITERATION

xPhys = x;
loop = 0;
change = 1;

  %% FE-ANALYSIS
 load vecs;load quad_d10_p5.mat;m=size(wdq,1);wquad=wdq;
 
 if nelx==4  vec=vec4; end;
 if nelx==10  vec=vec10; end;
 if nelx==20  vec=vec20; end;
 if nelx==50  vec=vec50; end;    
 
  cmean=0;dcmean=zeros(nely,nelx);
  cmean2=0;
  cdcmean=zeros(nely,nelx);
  svmmean=zeros(nely*nelx,1);
  sK = reshape(KE(:)*( [Emin+(xPhys(:)'.^penal).*(E0-Emin)]),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  
  F = sparse(2*(nely+1)*(nelx+1),1);
  ver_ind=fixeddofs(2:2:end)-(nely*2);
  hor_ind=fixeddofs(1:2:end)-(nely*2);
  F(ver_ind,1) = 2.0*0.01;
  F(hor_ind,1) = 0.0*0.01;
  U0 = zeros(2*(nely+1)*(nelx+1),1);
  U0(freedofs) = K(freedofs,freedofs)\F(freedofs);
  
  for i=1:10
      
      F = sparse(2*(nely+1)*(nelx+1),1);
      %ver_ind=fixeddofs(2:2:end)-(nely*2);
      hor_ind=fixeddofs(1:2:end)-(nely*2);
      %frc_hor = [ones(floor(length(hor_ind)/2),1); -ones(length(hor_ind)-floor(length(hor_ind)/2),1)];
      frc_hor=vec(:,i);
      %F(ver_ind,1) = 1.0;
      F(hor_ind,1) = frc_hor*0.01;
      U = zeros(2*(nely+1)*(nelx+1),1);
  %sK = reshape(KE(:)*( [Emin+(xPhys(:)'.^penal).*(E0-Emin)].*out(:,i)' ),64*nelx*nely,1);%%%%%%%%%%%%%
  
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  Umode(:,i)=U;
  end;
  UU=Umode*ndq'+U0*ones(1,m);
  for i=1:m
      U=UU(:,i);
%       ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);%%%%%%%%%%%%%%%%%%%%
%       c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
%       dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;dc=reshape(dc,nelx*nely,1);
      iter_cdc(:,i)=U;%[c;dc];      
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

