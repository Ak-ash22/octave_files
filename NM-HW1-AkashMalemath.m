clear all;
L = 0.328;        #length of the string
M_L = 0.660e-03;  #mass per length
F = 55.0;         #Force applied on the string
npt = 201;        #total number of space sample times
nst = 13;         #total number of time sample times
movie = 200;      #total number of frames
pinch = 2;
nmax = 20;        #max eigenvalue number
gamma = 100;      #damping factor


#Question 1
A = 1;
#For y = A*e^(-gamma*t)
#To get y(t_max) = y(t=0)/100
t_max = (1/100)*log(100);
fprintf('%6d\n',t_max);

#Question 2
xmin = 0;
xmax = L;
steps = (xmax-xmin)/(npt-1);
for i=1:npt
  x(i) = xmin + (i-1)*steps;
end

#Question 3
if pinch==1;
  p = L/2;
  z = 1;
  for i=1:npt
    if x(i)<=p;
      f(i) = z*x(i)/p;
    else;
      f(i) = z*(x(i)-L)/(p-L);
    endif
    g(i) = 0;
  endfor
elseif pinch==2;
  p1 = L/4;
  p2 = 3*L/4;
  z1 = 1;
  z2 = -1;
  for i=1:npt
    if x(i)<=p1;
      f(i) = z1*x(i)/p1;
    elseif x(i)>p1 && x(i)<=p2;
      f(i) = (z1-z2)*(x(i)-p1)/(p1-p2) + z1;
    else;
      f(i) = z2*(x(i)-L)/(p2-L);
    endif
    g(i) = 0;
  endfor
endif

##graph3 = plot(x,f,x,g);
##title('Graph for pinch = 1');
##xlabel('x');
##legend(f,g).show;
##waitfor(graph3);

#Question 4
for n=1:nmax
  for i=1:npt
    phi(n,i) = sqrt(2/L)*sin(n*pi*x(i)/L);
  endfor
end

##graph4 = plot(x,phi);
##title('nmax = 20 eigenfunctions of a vibrating string');
##xlabel('x');
##ylabel('phi(n,x)');
##waitfor(graph4);


#Question 5
c = sqrt(F/M_L);
for i=1:nmax
  n(i) = i;
  n_pi_L(i) = n(i)*pi/L;
  omega(i) = n_pi_L(i)*c;
  Omega(i) = sqrt(omega(i)^2 - gamma^2);
  v(i) = Omega(i)/(2*pi);
  T(i) = 1/v(i);
endfor

##fmt1 = '%4d %15.6E %15.6E %15.6E %15.6E\n';
##fprintf('   n        n_pi_L           Omega           v                T\n');
##fprintf('---------------------------------------------------------------------------\n');
##fprintf(fmt1,[n;n_pi_L;Omega;v;T]);

#Comparison of small omega and big omega for different values of damping factor gamma:
g1 = 10;
g2 = (10 + omega(1)/2)/2;
g3 = omega(1)/2;
gammas = [g1,g2,g3];
for i=1:length(gammas)
  #fprintf('   n     omega(n)          Omega(n)\n')
  for n=1:nmax
    n_pi_L(n) = n*pi/L;
    small_omega(n) = n_pi_L(n)*c;
    big_omega(n) = sqrt(small_omega(n)^2 - gammas(i)^2);
    #fprintf('%4d %15.6E %15.6E\n',n,small_omega(n),big_omega(n));
  endfor
endfor

#Question 6
for k=1:nmax
  I = phi(k,:).^2;
  ip(n) = trapz(x,I);
end
##fprintf('  n        <phi_n|phi_n>\n');
##fprintf('--------------------------\n');
##fprintf('%4d %15.3E\n',[n;ip]);


#Question 7
for k=1:nmax
  phi_f(k) = trapz(x,phi(k,:).*f);
  phi_g(k) = trapz(x,phi(k,:).*g);
end
##fprintf('    n      <phi_n|f>   <phi_n|g>\n');
##fprintf('----------------------------------\n');
##fprintf('%4d %15.3E %4d\n',[n;phi_f;phi_g]);


#Question 8
f_approx = 0;
g_approx = 0;
for k=1:nmax
  f_approx += phi_f(k)*phi(k,:);
  g_approx += phi_g(k)*phi(k,:);
end
##graph8 = plot(x,f_approx,x,g_approx);
##title('f_{approx}(x) and g_{approx}(x) vs x for nmax=20');
##xlabel('x');
##waitfor(graph8);


#Question 9
longest_period = max(T);
tmax = longest_period/2;
tmin = 0;
time_step = (tmax-tmin)/(nst-1);
for i=1:nst
  t(i) = tmin + (i-1)*time_step;
end

for i=1:nst
  for j=1:npt
    Psi(i,j) = 0;
    for n=1:nmax
      Psi(i,j) += phi(n,j)* (phi_f(n)*cos(Omega(n)*t(i)) + (1/Omega(n))*phi_g(n)*sin(Omega(n)*t(i)))*e^(-gamma*t(i));
    end
  end
##  graph9 = plot(x/L,Psi(i,:));
##  ts{k} = sprintf('t=%#10.2E s',t(i));
##  hold on;
end
##title('Amplitude patterns of Psi(x,t) for pinch=2');
##xlabel('x/L');
##ylabel('Psi(x,t)');
##set(lgd,'xlabel',' ','ylabel',' ');
##xlim([xmin/L,xmax/L]);
##ylim([-1,1]);
##waitfor(graph9);


#Question 10
normx = x/L;
z = 1;
tmin = 0;
tmax = 2*longest_period;
movie_time_step = (tmax-tmin)/(movie-1);
for i=1:npt
  psi(i) = 0;
  for n=1:nmax
    psi(i)+= phi(n,i) * phi_f(n);
  endfor
end

graph10 = plot(normx,psi,'XDataSource','normx','YDataSource','psi');
xlim([xmin/L,xmax/L]);
ylim([-z,z]);
line([xmin/L,xmax/L],[0,0],'linestyle','-','color',[0.5,0.5,0.5]);

for k=1:movie
  pause(0);
  t = tmin + movie_time_step*(k-1);
  for i=1:npt
    psi(i) = 0;
    for n=1:nmax
      psi(i) += phi(n,i)*(phi_f(n)*cos(Omega(n)*t) +(1/Omega(n))*phi_g(n)*sin(Omega(n)*t))*e^(-gamma*t);
    end
  end
  refreshdata();
end
waitfor(graph10);
