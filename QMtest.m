source('mystartdefaults.m');
tic

pos = get(groot,'DefaultFigurePosition');  #To get default matplotlib graphics

#Question 1
recipunit = 1.0E+10;
ekinscale = ((hbar*recipunit)^2/(2*elm))/qel;

#Question 2
xp_min = 0;
xp_max = 1;
n = 100;
step = (xp_max-xp_min)/n;


#Question 8
x_min = -10;
x_max = 10;
m = floor((x_max-x_min)/step);    #Calculating number of points between x_min and x_max with same step as inside perturbation


#Question 3
E_0 = 1;
k_0 = sqrt(E_0/ekinscale);
lambda = (2*pi/k_0);

xp = zeros(n,1);   #Making space discretization as a column
for i=1:n
  xp(i) = xp_min + step/2 + (i-1)*step;  #Discretization inside pertubation (Angstrom)
end

Phi0p = exp(1i*k_0*xp);  #Incident plane wave with amplitude 1 inside perturbation
#Phi0p is in column vector form


#Question 4
U0 = zeros(n,1);    #Initiating potential as a column vector
for i=1:n
  U0(i) = 2;  #Discretization of the perturbation
end
V0 = zeros(n,1);
for i=1:n
  V0(i) = U0(i)/ekinscale;  #Unit Conversion in k^2 units (1/Angstrom^2)
end


#Question 5
G0 = zeros(n,n);

for i=1:n
  for j=1:n
    G0(i,j) = step*exp(1i*k_0*abs(xp(i)-xp(j)))/(2i*k_0);
  endfor
end

#Question 6
T0 = eye(n,n) - G0*diag(V0);     #Scattering Matrix

#Question 7
PhiP = T0\Phi0p;

#check to verify PhiP
#hf(1) = figure('Question 7','off','name',['Square modulus of wavefunction inside perturbation'])
#g1 = plot(xp,abs(PhiP).^2);
#set(gcf,'Position',
##waitfor(g1);


#Question 8

##for i=1:20*n
##  x(i) = x_min + step/2 + (i-1)*step;
##end
##
##Phi0 = exp(1i*k_0*x');
##
##U = zeros(20*n,1);
##V = zeros(20*n,1);
##for i=1:20*n
##  U(i) = 2;
##  V(i) = U(i)/ekinscale;
##end
##
##G = zeros(20*n,20*n);
##for i=1:20*n
##  for j=1:20*n
##    G(i,j) = step*exp(1i*k_0*abs(x(i)-x(j)))/(2i*k_0);
##  endfor
##endfor
##
##
##T = eye(20*n,20*n) - G*diag(V);
##Phi = T\Phi0;
##
##g2 = plot(x,abs(Phi).^2);
##waitfor(g2);

#Question 9
x = zeros(m,1);
Phis = zeros(m,1);  #Phi scattered
Phi = zeros(m,1);
Prob = zeros(m,1);

for i=1:m
  x(i) = x_min + step/2 + (i-1)*step;
  Phis(i) = 0;
  for j=1:n
    Phis(i) += step*(exp(1i*k_0*abs(x(i)-xp(j))))/(2*1i*k_0) * V0(j) * PhiP(j);
  endfor
  Phi(i) = exp(1i*k_0*x(i)) + Phis(i);
  Prob(i) = abs(Phi(i))^2;
end
Ref = abs(Phis(1)/exp(1i*k_0*x(i)))^2   #Reflection  coeff
Tra = Prob(m)                            #Transmission Coeff
RpT = Ref+Tra                           #R + T = 1


#Question 10
w = (hbar*k_0^2)/(2*elm);
Psi = Phi*exp(1i*w*0);
speed_factor = recipunit*hbar/elm;  #in m/s
speed_factor = speed_factor*1e-6;   #in Mm/s
del_Psi = zeros(m-1,1);
for i=1:m-1
  del_Psi(i) = (Psi(i+1)-Psi(i))/(x(i+1)-x(i));
end
del_Psi(m) = del_Psi(m-1);
for i=1:m
  J(i) = speed_factor*real(conj(Psi(i))*del_Psi(i)/1i);
end

g3 = plot(x,abs(J));
waitfor(g3);


#Question 11
#Spectrum
Ugraph = zeros(m,1);
for i=1:m
  if (x(i)>=xp(1) && x(i)<=xp(n));
    for j=1:n
      if abs(x(i)-xp(j))<1e-12;
        Ugraph(i) = U0(j);
      endif
    endfor
  endif
end

hf(2) = figure('NumberTitle','off','name',['Scattering eigenstate']);
plot(x,Prob); hold on;
plot(x,real(Phi)); hold on;
plot(x, imag(Phi)); hold on;
plot(x,Ugraph); xlabel('x [Angstroms]');
legendstring{1} = '|\Phi(x)|^2';
legendstring{2} = 'Re(\Phi(x))';
legendstring{3} = 'Im(\Phi(x))';
legendstring{4} = 'U(x) [eV]';
lgs = legend(legendstring,'location','eastoutside');
##

#Question 12
##x_12_min = 0;
##x_12_max = 10;
##l = (x_12_max - x_12_min)/step;
##
##x_12 = zeros(l,1);
##phis = zeros(l,1);
##phi = zeros(l,1);
##prob = zeros(l,1);
##R = 0;
##T = 0;
##
##for i=1:m
##  x_12(i) = x_12_min + step/2 + (i-1)*step;
##  phis(i) = 0;
##  for j=1:n
##    phis(i) += step*exp(1i*k_0*abs(x_12(i)-xp(j)))/(2i*k_0) * V0(j) * PhiP(j);
##  endfor
##  phi(i) = exp(1i*k_0*x_12(i)) + phis(i);
##  prob(i) = abs(phi(i))^2;
##  R(i) = abs(phis(1)/(exp(1i*k_0*x_12(i))))^2;
##  T(i) = 1-R(i);
##end
##
##graph_12 = plot(x_12,R,x_12,T);
##waitfor(graph_12);

#Question 13
##E_13_min = 0;
##E_13_max = 10;
##l = 20;
##
##E_13 = zeros(l,1);
##phis = zeros(l,1);
##phi = zeros(l,1);
##prob = zeros(l,1);
##R = 0;
##T = 0;
##k = 0;
##
##for i=1:l
##
##  E_13(i) = E_13_min + (i-1)*step;
##  #phis(i) = 0;
##  k(i) = sqrt(E_13(i)/ekinscale);
##  for j=1:m
##    phis(j) = 0;
##    for s = 1:n
##      phis(j) += step*exp(1i*k(i)*abs(x(j)-xp(s)))/(2i*k(i)) * V0(s) * PhiP(s);
##    end
##    phi(j) = exp(1i*k(i)*x(i)) + phis(j);
##    prob(j) = abs(phi(j))^2;
##    R(i) = abs(phis(1)/(exp(1i*k(i)*x(i))))^2;
##    T(i) = 1-R(i);
##  endfor
##end
##
##graph_12 = plot(E_13,R,E_13,T);
##waitfor(graph_12);
toc
