A=load('A.dat');
B=load('B.dat');
C=load('C.dat');
F=load('F.dat');
X=load('X.dat');
n=length(F);
M=zeros(n,n);

for i=1:n/4
  ip1=mod(i,n/4)+1;
  im1=mod(i-2,n/4)+1;
  bloc=(i-1)*4+1;
  aloc=(im1-1)*4+1;
  cloc=(ip1-1)*4+1;
  M(bloc:bloc+3,bloc:bloc+3)=B(bloc:bloc+3,:);
  M(bloc:bloc+3,aloc:aloc+3)=A(bloc:bloc+3,:);
  M(bloc:bloc+3,cloc:cloc+3)=C(bloc:bloc+3,:);
end
Y=inv(M)*F;
X-Y

