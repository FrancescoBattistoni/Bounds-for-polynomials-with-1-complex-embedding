/* x = 0, t = 1 */

  {resz=
  y^11 - 1/3*y^10 + 43/81*y^9 - 43/243*y^8 - 4874/6561*y^7 + 4874/19683*y^6 - 
      4610/19683*y^5 + 4610/59049*y^4 + 2795/19683*y^3 - 2795/59049*y^2 - 
      575/19683*y + 575/59049
  ;}
  
  {res1g=
  y^4*z^3 - y^3*z^2 - 20/27*y^2*z^5 - 2/9*y^2*z^3 + 5/27*y^2*z + 16/27*y*z^4 + 
      10/27*y*z^2 + 1/27*y - 2/27*z^3 - 4/27*z
  ;}
  
  {Ly=
  y^2*z + 4/3*y*z*g - 2/3*y + 1/3*z - 2/3*g
  ;}
  
  R(y,z,g)=(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2);
  
  v=polrootsreal(resz,[0,1]);  /* Length 2 */
  y0=v[1];
  w=polrootsreal(subst(res1g,y,y0),[0,1]); /* Length 1 */
  z0=w[1];
  u=polrootsreal(subst(subst(Ly,y,y0),z,z0),[0,1]);  /* No roots */
  
  y0=v[2];
  w=polrootsreal(subst(res1g,y,y0),[0,1]); /* Length 3 */
  z0=w[1];
  u=polrootsreal(subst(subst(Ly,y,y0),z,z0),[0,1]);  /* No roots */
  
  z0=w[2];
  u=polrootsreal(subst(subst(Ly,y,y0),z,z0),[1/2,1]);  /* no roots */
  /* no stationary points */
....................................................
/* x = 0, t = 1, g = 1/2 */

L=(1+y+y^2)*(1-y*z)*(1+z+z^2)*(1-z+z^2)*sqrt(3)  <= 5.2
PROOF
*/ BORDER /*
* At y = 0
  (1+z+z^2)*(1-z+z^2)*sqrt(3) <= 3*sqrt(3) <= 5.2

* At z = 0
  (1+y+y^2)*sqrt(3) <= 3*sqrt(3) <= 5.2  

* At y = 1
  3*(1-z)*(1+z+z^2)*(1-z+z^2)*sqrt(3) <= 3*sqrt(3) <= 5.2  

* At z = 1
  (1+y+y^2)*(1-y)*3*sqrt(3) <= 3*sqrt(3) <= 5.2 

*/ INNER /*
L=(1+y+y^2)*(1-y*z)*(1+z+z^2)*(1-z+z^2)   /* times sqrt(3) */;
T=factor(deriv(L,y))
Ly=T[3,1];
T=factor(deriv(L,z))
Lz=T[1,1];
T=factor(polresultant(Ly,Lz,z))
res=T[2,1];

v=polrootsreal(res,[0,1]);     /* one root */
y0=v[1];
u=polrootsreal(subst(Ly,y,y0),[0,1]);     /* one root */
z0=u[1];

my(y=y0,z=z0); eval(L)*sqrt(3)  
\\ 3.2812203897729684509724761250679843614     <----------
END PROOF
....................................................

/* x = 1, t = 0 */

  {resz=
  y^11 - 1/2*y^10 + 421/375*y^9 - 421/750*y^8 + 542/1125*y^7 - 271/1125*y^6 + 
      106/5625*y^5 - 53/5625*y^4 - 287/5625*y^3 + 287/11250*y^2 - 23/1875*y + 
      23/3750
  ;}
  
  {res1g=
  y^5*z^2 - 4/5*y^4*z - 27/20*y^3*z^4 + 3/10*y^3*z^2 + 1/10*y^3 + 27/20*y^2*z^3 - 
      1/2*y^2*z - 1/4*y*z^2 + 1/5*y - 1/20*z
  ;}
  
  {Lz=
  y*z^2 + 4/3*y*z*g + 1/3*y - 2/3*z - 2/3*g
  ;}
  
  R(y,z,g)=(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*2*sqrt(1-g^2);
  
  v=polrootsreal(resz,[0,1]);  /* Length 2 */
  y0=v[1];
  w=polrootsreal(subst(res1g,y,y0),[0,1]); /* Length 2 */
  z0=w[1];
  u=polrootsreal(subst(subst(Lz,y,y0),z,z0),[0,1]);  /* No roots */
  z0=w[2];
  u=polrootsreal(subst(subst(Lz,y,y0),z,z0),[1/2,1]);  /* No roots */
  
  y0=v[2];
  w=polrootsreal(subst(res1g,y,y0),[0,1]); /* Length 1 */
  z0=w[1];
  u=polrootsreal(subst(subst(Lz,y,y0),z,z0),[1/2,1]);  /* no roots */
....................................................

/* x = 1, t = 0, g = 1/2 */

L = (1-y+y^2)*(1+y+y^2)*(1-y*z)*(1+z+z^2)*sqrt(3) <= 5.2
PROOF
*/ BORDER /*
* At y = 0
  (1+z+z^2)*sqrt(3) <= 3*sqrt(3) <= 5.2

* At z = 0
  (1-y+y^2)*(1+y+y^2)*sqrt(3) <= 3*sqrt(3) <= 5.2  

* At y = 1
  3*(1-z)*(1+z+z^2)*sqrt(3) <= 3*sqrt(3) <= 5.2  

* At z = 1
  (1-y+y^2)*(1+y+y^2)*(1-y)*3*sqrt(3) <= 3*sqrt(3) <= 5.2 

*/ INNER /*
L=(1-y+y^2)*(1+y+y^2)*(1-y*z)*(1+z+z^2) /* times sqrt(3) */;
T=factor(deriv(L,y))
Ly=T[2,1];
T=factor(deriv(L,z))
Lz=T[1,1];
T=factor(polresultant(Ly,Lz,z))
res=T[2,1];

v=polrootsreal(res,[0,1]);     /* one root */
y0=v[1];
u=polrootsreal(subst(Ly,y,y0),[0,1]);     /* one root */
z0=u[1];

my(y=y0,z=z0); eval(L)*sqrt(3)  
\\ 3.2812203897729684509724761250679843614     <----------
END PROOF
....................................................

