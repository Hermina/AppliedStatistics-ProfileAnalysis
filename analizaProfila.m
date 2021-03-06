gr1=load('grupa1.txt');
gr2=load('grupa2.txt');
gr3=load('grupa3.txt');
gr4=load('grupa4.txt');

sve=load('sve.txt');

B(1,:)=mean(gr1,1);
B(2,:)=mean(gr2,1);
B(3,:)=mean(gr3,1);
B(4,:)=mean(gr4,1);

plot(B')

lambda4(1,:)=ones(1,3);
lambda4(2:4,:)=-eye(3);

lambda11(1,:)=ones(1,10);
lambda11(2:11,:)=-eye(10);

J11=ones(11,1);
J4=ones(4,1);

novaGr1=gr1*lambda11;
novaGr2=gr2*lambda11;
novaGr3=gr3*lambda11;
novaGr4=gr4*lambda11;

s1=cov(gr1);
s2=cov(gr2);
s3=cov(gr3);
s4=cov(gr4);

e=5*s1+13*s2+14*s3+9*s4;
s=e/41;
ss=lambda11'*s*lambda11;
E1=lambda11'*e*lambda11

cv1=cov(novaGr1);
cv2=cov(novaGr2);
cv3=cov(novaGr3);
cv4=cov(novaGr4);
en=5*cv1+13*cv2+14*cv3+9*cv4;
sn=en/27;
X=diag([6,14,15,10]);

H=(lambda4'*B*lambda11)'*inv(lambda4'*inv(X)*lambda4)'*lambda4'*B*lambda11
hip=det(E1+H);
pog=det(E1);
rez=pog/hip

v = 10*(4-1);
X2 = (-1)*(((44)-(.5*(10+4)))*log(rez));
P = 1-chi2cdf(X2,v)


H2=(lambda4'*B*J11)'*inv(lambda4'*inv(X)*lambda4)'*lambda4'*B*J11;
E2=J11'*e*J11;
hip2=det(E2+H2);
pog2=det(E2);
rez2=pog2/hip2
v = 3;
X22 = (-1)*(((44)-(.5*(10+4)))*log(rez2));
P2 = 1-chi2cdf(X22,v)

H3=(J4'*B*lambda11)'*inv(J4'*inv(X)*J4)'*J4'*B*lambda11;
E3=lambda11'*e*lambda11;
hip3=det(E3+H3);
pog3=det(E3);
rez3=pog3/hip3

v = 10*(4-1);
X23 = (-1)*(((44)-(.5*(10+4)))*log(rez3));
P3 = 1-chi2cdf(X23,v)

