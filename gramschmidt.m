a = [1,4,6;
    -1,7,2;
     4,1,9];

y = 1e-10;

b = [1,1,1;
     y 0 0;
     0 y 0;
     0 0 y];

c = hilb(10);

x=b;
[s,t]=GramSchmidtEuclidean(x,'stable',false,'normalized',true);
[q,r]=GramSchmidtEuclidean(x,'stable',true,'normalized',true);

m=rows(q); n = columns(q);
if (n<=10)
  printtable(x,'precision',3,'Title','Input matrix')
  printtable(s,'precision',3,'Title','q from classical algorithm (unstable)')
  printtable(t,'precision',3,'Title','r from classical algorithm (unstable)')
  printtable(q,'precision',3,'Title','q from modified algorithm (stable)')
  printtable(r,'precision',3,'Title','r from modified algorithm (stable)')
  printtable(s'*s,'precision',3,'Title','q''*q from classical algorithm (unstable)')
  printtable(q'*q,'precision',3,'Title','q''*q from modified algorithm (stable)')
  printtable(s*t,'precision',3,'Title','q*r from classical algorithm (unstable)')
  printtable(q*r,'precision',3,'Title','q*r from modified algorithm (stable)')
  printtable(s*t-x,'precision',3,'Title','q*r-a from classical algorithm (unstable)')
  printtable(q*r-x,'precision',3,'Title','q*r-a from modified algorithm (stable)')
end

result1 = norm(s*t-x);
result2 = norm(q*r-x);

disp(['|q*r-a| for classical algorithm (unstable) is ' num2str(result1)]);
disp(['|q*r-a| for modified algorithm (stable) is ' num2str(result2)]);

result3 = norm(diag(s'*s)-s'*s);
result4 = norm(diag(q'*q)-q'*q);
disp(['|diag(q''*q)-q''q| for classical algorithm (unstable) is ' num2str(result3)]);
disp(['|diag(q''*q)-q''q| for modified algorithm (stable) is ' num2str(result4)]);
