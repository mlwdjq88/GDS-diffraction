function Z = ZgenNM(n,m)
R = RgenNM(n,abs(m));
A = AgenNM(m);
Z = @(r, phi) R(r).*A(phi);


function f = RgenNM(n,m)
p = (n-m)/2;
f = @(r) 0;

for k = 0:p
    coef = (-1)^k*nchoosek(n-k,k)*nchoosek(n-2*k,(n-m)/2-k);
    ex = (n - 2*k);
    f = @(r) f(r) + coef*r.^ex;
end


% Forms azimuthal function based on n,m
function g = AgenNM(m)
if m > 0
    g = @(phi) cos(m*phi);
elseif m < 0
    g = @(phi) -sin(m*phi);
else
    g = @(phi) 1;
end