
function dy = eq_nondim_K(t,y)



global rn ri rm alphani alphanm alphamn alphain eps Kn Ki Km epsi

n = y(1);
i0 = y(2);
im = y(3);
dy = zeros(3,1);

dy(1) = rn*n  - rn*n^2/Kn         - alphain*i0*n   - alphamn*im*n;
dy(2) = ri*i0 - ri*i0*(im+i0)/Ki  - alphani*n*i0   - eps*n*i0 - epsi*im*i0;
dy(3) = rm*im - rm*im*(im+i0)/Km  - alphanm*n*im   + eps*i0*n + epsi*i0*im;

end