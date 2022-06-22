a = zeros(1,5);

ntemp = 6;
imt = 5;
ai = zeros(1,ntemp);
for i = 1:imt-1
    ai(i) = a(i);
end
ai(imt) = 1;
ai(imt + 1) = 2;

for i = imt+2:ntemp
    ai(i) = a(i-1);
end
    
