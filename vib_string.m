#Question 1
clear all;
L = 0.328;
mass_per_length = 0.6603e-3;
F = 55.0;
npt = 201;
nst = 13;
pinch_defined = 2;

#Question 2
xmin = 0;
xmax = L;
step = (xmax-xmin)/(npt-1);
for i=1:npt
  x(i) = xmin + (i-1)*step;
end

#Question 3

if pinch_defined==0;
  printf('No pinch given\n')
  f = [0]*length(x);
  g = 0;
elseif pinch_defined==1;
  p = L/2;
  z = 1;
  for i=1:npt
    if x(i)<p;
      f(i) = z*x(i)/p;
      g(i) =0;
    else;
      f(i) = z*(x(i)-L)/(p-L);
      g(i) =0;
    endif
  endfor
elseif pinch_defined==2;
  p1 = L/4;
  p2 = 3*L/4;
  z1 = 1;
  z2 = -1;
  for i=1:npt
    if x(i)<=p1;
      f(i) = z1*x(i)/p1;
      g(i)=0;
    elseif x(i)>p1 && x(i)<=p2;
      f(i) = (z1-z2)*(x(i)-p1)/(p1-p2) + z1;
      g(i)=0;
    else;
      f(i) = z2*(x(i)-L)/(p2-L);
      g(i)=0;
    endif
  endfor
endif

#pinch_graph = plot(x,f);
#waitfor(pinch_graph);


#Question 4
nmax = 20;
for k=1:nmax
  for i=1:npt
    eig_fun(k,i) = sqrt(2/L)*sin(k*pi*x(i)/L);
  end
  #eig_graph = plot(x,eig_fun);
  #hold on
end

#Question 5
c = sqrt(F/mass_per_length);

for i=1:nmax
  n(i) = i;
  n_pi_L(i)  = n(i)*pi/L;
  w(i) = n_pi_L(i)*c;
  v(i) = w(i)/(2*pi);
  T(i) = 1/v(i);
endfor

##fmt = '%15.1E %15.3E %15.3E %15.3E %15.3E\n';
##fprintf('           n             n_pi_L           w               v               T\n')
##fprintf('   -------------------------------------------------------------------------------\n')
##fprintf(fmt,[n;n_pi_L;w;v;T])

#T = table(n',n_pi_L',w',v',T','VariableNames',{'n','n*pi/L','w','v','T'})

#Question 6

for i=1:nmax
  I = eig_fun(i,:).^2;
  inner_prod(i) = trapz(x,I);
end

fmt1 = '%15.1E %15.3E\n';
fprintf('        n            <phi(n)|phi(n)>\n')
fprintf('     ----------------------------------\n')
fprintf(fmt1,[n;inner_prod])


#Question 7

for k=1:nmax
  phi_f(k) = trapz(x,eig_fun(k,:).*f);
  phi_g(k) = trapz(x,eig_fun(k,:).*g);
end
fprintf('    n      <phi_n|f>   <phi_n|g>\n');
fprintf('----------------------------------\n');
fprintf('%4d %15.3E %4d\n',[n;phi_f;phi_g]);


#Question 8
f_approx = 0;
g_approx = 0;
for n=1:nmax
  f_approx = f_approx + phi_f(n)*eig_fun(n,:);
  g_approx = g_approx + phi_g(n)*eig_fun(n,:);
endfor
##plot(x,f_approx,x,g_approx);
##hold on;

#Question 9
longest_period = max(T);
nmax=20;
nst = 7;
movie = 200;
tmin = 0;
tmax = longest_period/2;
time_step = (tmax-tmin)/(nst-1);
for i=1:nst
  t(i) = tmin + (i-1)*time_step;
endfor


for j=1:nst
  for i=1:npt
    Psi(j,i)=0;
    for n=1:nmax
      Psi(j,i) = Psi(j,i) + eig_fun(n,i)*phi_f(n)*cos(w(n)*t(j));
    endfor
  endfor
##  graph9 = plot(x/L,Psi(j,:));
##  ts{k} = sprintf('t=%#10.2E s',t(j));
##  hold on;
endfor
#lgd = legend(ts,'location','eastoutside');
###set(lgd,'xlabel',' ','ylabel',' ');
##xlim([xmin/L,xmax/L]);
##ylim([-1,1]);
##waitfor(graph9);





#Question 10
normx = x/L;
z = 1;
tmax = longest_period*2;
time_step1 = (tmax-tmin)/(movie-1);
for i=1:npt
  psi(i)=0;
  for n=1:nmax
    psi(i)+= eig_fun(n,i)*phi_f(n);
  endfor
end
g = plot(normx,psi,'XDataSource','normx','YDataSource','psi');
xlim([0,1]);
ylim([-z,z]);
line([0,1],[0,0],'linestyle','-','color',[0.5,0.5,0.5]);

for k=1:movie
  pause(0);
  t = tmin +(k-1)*time_step1;
  for i=1:npt
    psi(i)=0;
    for n=1:nmax
      psi(i) += eig_fun(n,i)*(phi_f(n)*cos(w(n)*t) + (1/w(n))*phi_g(n)*sin(w(n)*t))*e^(-100*t);
    endfor
  endfor
  refreshdata();
endfor
waitfor(g);

