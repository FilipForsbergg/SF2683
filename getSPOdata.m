
function [lambdavec,Tvec,cvec] = getSPOdata(birthdate)   %yyyymmdd
    
lambdavec = 1/1000*[50 40 45 50 25 48 60 35 15];
Tvec      = [6 8 14 25 12 18 33 8 12];
cvec      = [12 14 21 20 11 45 75 30 22];
%Cmax      = 500;
n=length(lambdavec);

rng(birthdate)
lambdavec = lambdavec.*(1+0.2*rand(1,n)); 
Tvec      = round(Tvec.*(1+0.2*rand(1,n))); 
cvec      = round(cvec.*(1+0.2*rand(1,n))); 

end