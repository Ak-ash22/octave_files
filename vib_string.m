#Question 1
clear all;
L = 0.328;
mass_per_length = 0.6603e-3;
F = 55.0;
npt = 201;
nst = 13;

#Question 2
xmin = 0;
xmax = L;
step = (xmax-xmin)/(npt-1);
for i=1:npt
  x(i) = xmin + (i-1)*step;
end

#Question 3
pinch_defined = 1;

for k=1:length(pinch_defined)
  if pinch_defined(k)==0;
    printf('No pinch given\n')
    f0 = [0]*length(x);
  endif

  if pinch_defined(k)==1;
    p = L/2;
    z = 1;
    for i=1:npt
      if x(i)<p;
        f1(i) = z*x(i)/p;
      else;
        f1(i) = z*(x(i)-L)/(p-L);
      endif
    endfor
  endif

  if pinch_defined(k)==2;
    p1 = L/4;
    p2 = 3*L/4;
    z1 = 1;
    z2 = -1;
    for i=1:npt
      if x(i)<=p1;
        f2(i) = z1*x(i)/p1;
      elseif x(i)>p1 && x(i)<=p2;
        f2(i) = (z1-z2)*(x(i)-p1)/(p1-p2) + z1;
      else;
        f2(i) = z2*(x(i)-L)/(p2-L);
      endif
    endfor
  endif
endfor
#pinch_graph = plot(x,f0,x,f1,x,f2);
#waitfor(pinch_graph);


#Question 4
#eig_no = linspace(1,10,10);
for k=1:20
  for i=1:npt
    eig_fun(k,i) = sqrt(2/L)*sin(k*pi*x(i)/L);
  end
  #eig_graph = plot(x,eig_fun);
  #hold on
end

#Question 5
c = sqrt(F/mass_per_length);
n = linspace(1,20,20);

for i=1:20
  n_pi_L(i)  = n(i)*pi/L;
  w(i) = n_pi_L(i)*c;
  v(i) = w(i)/(2*pi);
  T(i) = 1/v(i);
endfor

fmt = '%15.1E %15.3E %15.3E %15.3E %15.3E\n';
fprintf('           n             n_pi_L           w               v               T\n')
fprintf('   -------------------------------------------------------------------------------\n')
fprintf(fmt,[n;n_pi_L;w;v;T])

#T = table(n',n_pi_L',w',v',T','VariableNames',{'n','n*pi/L','w','v','T'})

#Question 6
y = 0:L/100:L;
for i=1:20
  phi = sqrt(2/L)*sin(i*pi*y/L);
  I = phi.^2;
  inner_prod(i) = trapz(y,I);
end

#fmt1 = '%15.1E %15.3E\n';
#fprintf('        n            <phi(n)|phi(n)>\n')
#fprintf('     ----------------------------------\n')
#fprintf(fmt1,[eig_no;inner_prod])


#Question 7
pinch = [1,2];
for k=1:2
  for j=1:20
    if k==1;
      p = L/2;
      r1 = 0:p/100:p;
      r2 = p:(L-p)/100:L;
      z = 1;
      f_1 = z*r1/p;
      f_2 = z*(r2-L)/(p-L);
      phi_n1 = sqrt(2/L)*sin(j*pi*r1/L);
      phi_n2 = sqrt(2/L)*sin(j*pi*r2/L);
      I_1 = phi_n1.*f_1;
      I_2 = phi_n2.*f_2;
      ip_1(j) = trapz(r1,I_1);
      ip_2(j) = trapz(r2,I_2);
      ip1(j) = ip_1(j)+ip_2(j);


    else;
      p1 = L/4;
      p2 = 3*L/4;
      r1 = 0:p1/100:p1;
      r2 = p1:(p2-p1)/100:p2;
      r3 = p2:(L-p2)/100:L;
      z1 = 1;
      z2 = -1;
      f_1 = z1*r1/p1;
      f_2 = (z1-z2)*(r2-p1)/(p1-p2) + z1;
      f_3 = z2*(r3-L)/(p2-L);
      phi_n1 = sqrt(2/L)*sin(j*pi*r1/L);
      phi_n2 = sqrt(2/L)*sin(j*pi*r2/L);
      phi_n3 = sqrt(2/L)*sin(j*pi*r3/L);
      I_1 = phi_n1.*f_1;
      I_2 = phi_n2.*f_2;
      I_3 = phi_n3.*f_3;
      ip_1(j) = trapz(r1,I_1);
      ip_2(j) = trapz(r2,I_2);
      ip_3(j) = trapz(r3,I_3);
      ip2(j) = ip_1(j)+ip_2(j)+ip_3(j);
    endif
  endfor
endfor
#fmt2 = '%4d %15.6G\n';
#fprintf('For 1 Pinch\n');
#fprintf('        n             <phi_n|f>\n');
#fprintf('      ---------------------------\n');
#fprintf(fmt2,[eig_no;ip1]);
#fprintf('For 2 Pinches\n');
#fprintf('        n             <phi_n|f>\n');
#fprintf('      ---------------------------\n');
#fprintf(fmt2,[eig_no;ip2])


#Question 7

##f_approx = 0;
##for n=1:10
##  f_approx = f_approx .+ trapz(x,eig_fun(:,n).*f1)*eig_fun(:,n)
##endfor
##
##plot(x,f_approx);
##hold on;







#Question 8
f_approx = 0;
for n=1:20
  f_approx = f_approx + trapz(x,eig_fun(n,:).*f1) .*eig_fun(n,:);
endfor
#plot(x,f_approx);
#hold on;

#Question 9
longest_period = T(1);
nmax=20;
nst = 7;
movie = 200;
tmin = 0;
tmax = 2*longest_period;
time_step = (tmax-tmin)/(movie-1);
for i=1:nst
  t(i) = tmin + (i-1)*time_step;
endfor



for i=1:npt
  psi(i)=0;
  for n=1:nmax
    psi(i) = psi(i) + eig_fun(n,i).*trapz(x,eig_fun(n,:).*f1);
  endfor
endfor
plot(x,psi);
hold on;

#Question 10
normx = x/L;
z = 1;


#g = plot(normx,psi,'XDataSource','normx','YDataSource','psi');
xlim([0,1]);
ylim([-z,z]);
line([0,1],[0,0],'linestyle','-','color',[0.5,0.5,0.5]);

for k=1:movie
  pause(0);
  t = tmin +(k-1)*time_step;
  for i=1:npt
    psi(i)=0;
    for n=1:nmax
      psi(i) = psi(i) + eig_fun(n,i)*trapz(x,eig_fun(n,:).*f1)*cos(w(n)*t);
    endfor
  endfor
  refreshdata();
endfor
#waitfor(g);

