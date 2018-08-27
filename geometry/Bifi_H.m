function [trajectory_compliance,trajectory_design,trajectory_gradient]=Bifi_H(nelx,nely,volfrac,rmin,penal)

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

x = repmat(volfrac,nely,nelx);
x(1:nely/2,nelx/2:nelx)=0.001;
xPhys = x;
loop = 0;
change = 1;
%% START ITERATION
while ((change > 0.01) && (loop<277))
  loop = loop + 1;
  %% FE-ANALYSIS
  load data_100.mat;load des_quad.mat;wquad=w';  
  m=size(wquad,2);
  xfilt=(H*reshape(xPhys,nelx*nely,1))./Hs;

  
  %%%%%%%%%%%%%Bifi%%%%%%%%%%
   
   mesh_l=10;zrb=100/mesh_l;
   
    for i=1:mesh_l
        for j=1:mesh_l
            xL(i,j)=sum(sum(xPhys(1+(i-1)*zrb:i*zrb,1+(j-1)*zrb:j*zrb)))/(zrb^2);
        end;
    end;
    
   [UL]=randomize_L(mesh_l,mesh_l,3,3,mesh_l/2+1,xL);
   [a11,a22,perm]=qr(UL,'vector');
   %imp=rank(UL)
   %size(UL)
   imp=10;
   for i=1:(m-imp)
       unk=UL(:,perm(imp+i));
       coef_L(:,i)=(UL(:,perm(1:imp))'*UL(:,perm(1:imp)))\(UL(:,perm(1:imp))'*unk);
   end;
  
  [UH_imp]=randomize_H(100,100,3,6,50,xPhys,perm(1:imp));
  if imp==43
      UH=UH_imp;
  else
      UH=[UH_imp UH_imp*coef_L];
  end;
  
  for i=1:m
      indperm=find(perm==i);
      UHH(:,i)=UH(:,indperm);
  end;
  

  cmean=0;dcmean=zeros(nely,nelx);dvmean=dcmean;
  cmean2=0;
  cdcmean=zeros(nely,nelx);
  eta=out_100*0.1+0.45;beta=2;
  for i=1:m
      ind=perm(i);ind=i;
      U=UHH(:,ind);
      for j=1:(nelx*nely)
          xthr(j,1) = (tanh(beta*eta(j,ind)) + tanh(beta*(xfilt(j) - eta(j,ind))) ) /  (tanh(beta*eta(j,ind)) + tanh(beta*(1 - eta(j,ind))) );
          deri_x(j,1) = ( beta * (sech(beta*(xfilt(j)-eta(j,ind))))^2 ) / (tanh(beta*eta(j,ind)) + tanh(beta*(1 - eta(j,ind))) );
      end;
      xPhys=reshape(xthr,nely,nelx);
      
      %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
      ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);%%%%%%%%%%%%%%%%%%%%
      c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
      dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
      
      dc=(H'* [deri_x.*reshape(dc,nelx*nely,1)]   )./Hs;dc=reshape(dc,nely,nelx);
      dv=reshape((deri_x),nely,nelx);dvmean=dvmean+dv*wquad(i);
      
      
      
      cmean=c*wquad(ind)+cmean;
      cmean2=c^2*wquad(ind)+cmean2;
      
      dcmean=dcmean+dc*wquad(ind);
      cdcmean=cdcmean+c*dc*wquad(ind);
      
      
  end;
  c=cmean+0.1*sqrt(cmean2-cmean^2);dc=dcmean+0.1*([cdcmean-cmean*dcmean]/sqrt(cmean2-cmean^2));
  
  dv=dvmean;
 
  %% FILTERING/MODIFICATION OF SENSITIVITIES
%   if ft == 1
%     dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
%   elseif ft == 2
%     dc(:) = H*(dc(:)./Hs);
%     dv(:) = H*(dv(:)./Hs);
%   end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    xnew(1:nely/2,nelx/2:nelx) = 0.001;
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  
     
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
  mean(xPhys(:)),change);
  %% PLOT DENSITIES
  figure(1)
   colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;

trajectory_design(:,:,loop)=xPhys;
trajectory_gradient(:,:,loop)=dc;
trajectory_compliance(loop)=c;
  
end

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

