Temp = [400 800];




fsfit =  [0.0075   -0.1993    0.3305    0.0896    0.0004   -0.0005    0.5187    0.6120    0.0019   -0.1187    0.4251    0.1508 0.0038   -0.2082    0.2238    0.0812];

A = [   fsfit(1:4);
        fsfit(5:8);
        fsfit(9:12);
        fsfit(13:16)];
    
AA = sortrows(A,3);



ind = [3 7 11 15];
fsfit(ind) = fsfit(ind)*(Temp(end)-Temp(1)) + Temp(1); 
ind = [4 8 12 16];
fsfit(ind)  = fsfit(ind)*(Temp(end)-Temp(1));
