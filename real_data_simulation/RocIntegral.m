function value = RocIntegral(pd, pf, k)

value = pd(1)*pf(1);
for i = 1:k
value = value + (pd(i+1) + pd(i))*(pf(i+1) - pf(i))/2;
end
end
 
