function I = innerint2(t,T,E)
I = trapz(t,exp(-1000*E./(8.31446261815324*T)));
end