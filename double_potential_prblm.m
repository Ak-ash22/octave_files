source('mystartdefaults.m');
pos = get(groot,'DefaultFigurePosition');  #To get default matplotlib graphics
tic

#Declaring units
recipunit = 1.0E+10;
ekinscale = ((hbar*recipunit)^2/(2*elm))/qel;

#Energy E0
E0_min = 0;
E0_max = 0.3;
l = 5;
estep = 0.1;
E0 = 0;
for i=1:l
  E0(i) = 0.15;
end

#space Discretization inside potentials

#Inside U1
x1_min = 0;
x1_max = 80;
n = 400;

#Inside U2
step = (x1_max-x1_min)/n;

U = zeros(n,1);    #Potential barriers;

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
x_min = -20;
x_max = 100;
m = floor((x_max-x_min)/step);
x = zeros(m,1);
for i=1:m
  x(i) = x_min + step/2 + (i-1)*step;
end

#Defining parameters
k = 0;
for i=1:l
  k(i) = sqrt(E0(i)/ekinscale);
  lmbda(i) = 2*pi/k(i);
end

#plane wave inside potentials
phi1 = {};
for i=1:l
  phi1{i} = exp(1i*k(i)*x1);
end

#Greens function
G1 = {};
for i=1:l
  G1{i} = zeros(n,n);
  for j=1:n
    for h=1:n
      G1{i}(j,h) = step*exp(1i*k(i)*abs(x1(j)-x1(h)))/(2*1i*k(i));
    end
  end
end

#Scattering matrix and scattering eigenfunction inside perturbation
T1 = {};
Phip1 = {};
for i=1:l
  T1{i} = eye(n,n) - G1{i}*diag(V);
  Phip1{i} = T1{i}\phi1{i};
end

graph1= plot(x1,abs(Phip1{3}).^2);
waitfor(graph1);

#Working outside potential : To find phi scattered

Phis = {};
Phi = {};
Prob = {};
for i=1:l
  Phis{i} = zeros(m,1);
  Phi{i} = zeros(m,1);
  Prob{i} = zeros(m,1);
  for j=1:m
    Phis{i}(j) = 0;
    for h=1:n
       Phis{i}(j) += step*exp(1i*k(i)*abs(x(j)-x1(h)))/(2*1i*k(i))*V(h)*Phip1{i}(h);
    endfor
    Phi{i}(j) = exp(1i*k(i)*x(j)) + Phis{i}(j);
    Prob{i}(j) = abs(Phi{i}(j))^2;
    #R{i}(j) = abs(Phis{i}(j)/exp(1i*k(i)*x(j)))^2;
  endfor

  #T{i}(j) = Prob{i}(m);
  #total{i} = R{i} + T{i};
end


#Spectrum
Ugraph = zeros(m,1);
for i=1:m
  if (x(i)>=x1_min && x(i)<=15);
    Ugraph(i) = 0.2;
  elseif (x(i)>65 && x(i)<=x1_max);
    Ugraph(i) = 0.1;
  endif
end

hf(2) = figure('NumberTitle','off','name',['Scattering eigenstate']);
plot(x,Prob{3}); hold on;
plot(x,real(Phi{3})); hold on;
plot(x, imag(Phi{1})); hold on;
plot(x,Ugraph); xlabel('x [Angstroms]');
legendstring{1} = '|\Phi(x)|^2';
legendstring{2} = 'Re(\Phi(x))';
legendstring{3} = 'Im(\Phi(x))';
legendstring{4} = 'U(x) [eV]';
lgs = legend(legendstring,'location','eastoutside');

