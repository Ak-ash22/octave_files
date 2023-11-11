source('mystartdefaults.m');
pos = get(groot,'DefaultFigurePosition');  #To get default matplotlib graphics
tic

#Declaring units
recipunit = 1.0E+10;
ekinscale = ((hbar*recipunit)^2/(2*elm))/qel;

#Energy E0
E0 = 0:0.005:0.3;
Ref = zeros(length(E0),1);
Tra = zeros(length(E0),1);
RpT = zeros(length(E0),1);
#space Discretization inside potentials

#Inside U1
x1_min = 0;
x1_max = 80;
n = 80;

#Inside U2
step = (x1_max-x1_min)/n;

U = zeros(n,1);    #Potential barrier;
x1 = zeros(n,1);

for i=1:n
  x1(i) = x1_min + step/2 + (i-1)*step;
  if x1(i)<=15
    U(i) = 0.2;
  elseif x1(i)>=65 && x1(i)<=80
    U(i) = 0.1;
  else
    U(i) = 0;
  endif
  V(i) = U(i)/ekinscale;
end

#Space discretization
##x_min = -20;
##x_max = 100;
##m = floor((x_max-x_min)/step);
##x = zeros(m,1);
##for i=1:m
##  x(i) = x_min + step/2 + (i-1)*step;
##end

#Defining parameters
k_0 = 0;
for i=1:length(E0)
  k_0(i) = sqrt(E0(i)/ekinscale);
  lmbda(i) = 2*pi/k_0(i);
end

m=2;
x = zeros(m,1);
Phis = zeros(m,1);
Phi = zeros(m,1);
x(1) = x1_min-1;
x(2) = x1_max+1;

fprintf('  k             E0                R             T              R+T\n');
for k=1:length(E0)
    if E0(k)==0
      Ref(k)=1;
      Tra(k)=0;
      RpT(k)=1;
    else
      Phi0p = exp(1i*k_0(k)*x1);
      G0 = zeros(n,n);
      for j=1:n
        for i=1:n
          G0(i,j) = step*exp(1i*k_0(k)*abs(x1(i)-x1(j)))/(2i*k_0(k));
        endfor
      endfor

      T = eye(n,n)-G0*diag(V);
      Phip = T\Phi0p;

      for i=1:m
        Phis(i) = 0;
        for j=1:n
          Phis(i) += step*(exp(1i*k_0(k)*abs(x(i)-x1(j))))/(2*1i*k_0(k)) * V(j) * Phip(j);
        endfor
        Phi(i) = exp(1i*k_0(k)*x(i)) + Phis(i);
      endfor
      Ref(k) = abs(Phis(1)/exp(1i*k_0(k)*x(1)))^2;
      Tra(k) = abs(Phi(2))^2;
      RpT(k) = Ref(k)+Tra(k);
    endif
    #fprintf('%15.3E %15.3E %15.3E %15.3E %15.3E\n',k_0(k),E0(k),Ref(k),Tra(k),RpT(k));
endfor

spec_graph = plot(E0,Ref,E0,Tra);
waitfor(spec_graph);
