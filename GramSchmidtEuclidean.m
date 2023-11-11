#Gram Schmidt
function[q,r] = GramSchmidtEuclidean(a,varargin)

defaultStable = true;
defaultNormalized = true;

p = inputParser; p.CaseSensitive = false; %initializee parser

%Input validation functions
validMatrix = @(x) ismatrix(x) && isreal(x);
validbool = @(x) islogical(x);

%Define rules
addRequired(p,'a',validMatrix);
addParameter(p,'normalized',defaultNormalized,validbool);
addParameter(p,'stable',defaultStable,validbool);

parse(p,a,varargin{:});
a = p.Results.a;
stable = p.Results.stable;
normalized = p.Results.normalized;

[n,k] = size(a);
q = zeros(n,k);
r = zeros(k,k);

if (n<k)
  error(['Number of rows of input matrix must be > than its number of columns'])
end

if (~stable && ~normalized)
  for j=1:k
    q(:,j) = a(:,j);
    for i=1:j-1
      r(i,j) = q(:,i)'*a(:,j)/(q(:,i)'*q(:,i));
      q(:,j) = a(:,j) - r(i,j)*q(:,j);
    endfor
    r(j,j) = q(:,j)'*q(:,j)/(q(:,j)'*q(:,j));
  endfor
end


if (~stable && normalized)
  for j=1:k
    q(:,j) = a(:,j);
    for i=1:j-1
      r(i,j) = q(:,i)'*a(:,j);
      q(:,j) = a(:,j) - r(i,j)*q(:,j);
    endfor
    r(j,j) = q(:,j)'*q(:,j);
  endfor
end

if (stable && ~normalized)
  for j=1:k
    q(:,j) = a(:,j);
  end
  for i=1:k
    r(i,i)=1;
    for j=i+1:k
      r(i,j) = q(:,i)'*q(:,j)/(q(:,i)'*q(:,i));
      q(:,j) = q(:,j) - r(i,j)*q(:,i);
    endfor
  endfor
end

if (stable && normalized)
  for i=1:k
    q(:,i) = a(:,i);
  endfor
  for i=1:k
    r(i,i) = norm(q(:,i));
    q(:,i) = q(:,i)/r(i,i);
    for j=i+1:k
      r(i,j) = q(:,i)'*q(:,j);
      q(:,j) = q(:,j) - r(i,j)*q(:,i);
    endfor
  endfor
end

