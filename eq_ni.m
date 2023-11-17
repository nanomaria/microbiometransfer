function dy = eq_ni(t,y)

global Kn Ki rn ri alphain alphani


n = y(1);
i0 = y(2);
dy = zeros(2,1);


dy(1) = rn* n * (1-n/Kn)  - alphain*n*i0;
dy(2) = ri* i0* (1-i0/Ki) - alphani*i0*n; 


end