
/* x = 1 */
  L=2*(1-y^2)*(1-(y*z)^4)*(1+(y*z)^2)

  Ty=factor(deriv(L,y));
  Tz=factor(deriv(L,z));
  Ly=Ty[3,1];
  Lz=Tz[6,1];

  T=factor(polresultant(Ly,Lz,y))
  /* no roots */
  /* L <= 2*/
.........................................................................

/* z = 1, t = 0 */
  L=(1+x)*(1-(x*y)^2)*(1-y^2)

  Tx=factor(deriv(L,x));
  Ty=factor(deriv(L,y));
  Lx=Tx[3,1];
  Ly=Ty[3,1];

  T=factor(polresultant(Lx,Ly,x))
  /* no roots */
  /* L <= 2*/
.........................................................................

/* z = 1, t = 1  */
  L=(1+x)*(1+x*y)*(1-x*y)*(1-2*g*x*y+(x*y)^2)*(1-y)*(1+y)*(1+2*g*y+y^2)
  
  Tx=factor(deriv(L,x));
  Ty=factor(deriv(L,y));
  Tg=factor(deriv(L,g));
  Lx=Tx[4,1];
  Ly=Ty[2,1];
  Lg=Tg[7,1];
  
  T=factor(polresultant(Lx,Lg,g))
  res1g=T[3,1]*T[4,1]
  T=factor(polresultant(Ly,Lg,g))
  res2g=T[5,1]
  
  T=factor(polresultant(res1g,res2g,y))
  resy=T[3,1]
  
  x0=5/9;
  
  polrootsreal(subst(res1g,x,x0),[0,1])
  y0 = 0.92338051687663868591100756630627792309
  
  polrootsreal(subst(subst(Lg,x,x0),y,y0),[0,1])
  g0 = 0.24439036756006639466066658614797554712
  
  my(x=x0,y=y0,g=g0);(1+x)*(1+x*y)*(1-x*y)*(1-2*g*x*y+(x*y)^2)*(1-y)*(1+y)*(1+2*g*y+y^2)
  \\ 0.39400332150086052998582522674469259849    <---- value
  
/* g = 0
    (1+x)*(1-(x*y)^4)*(1-y^4)  <= 2 */

/* g = 1
  (1+x)*(1+x*y)*(1-x*y)^3*(1-y)*(1+y)^3  */
  L=(1+x)*(1+x*y)*(1-x*y)^3*(1-y)*(1+y)^3;
  Tx=factor(deriv(L,x));
  Ty=factor(deriv(L,y));
  Lx=Tx[4,1];
  Ly=Ty[4,1];

  T=factor(polresultant(Lx,Ly,y))
  res1g=T[3,1]*T[4,1]

  T=factor(polresultant(res1g,res2g,y))
  resy=T[3,1]

  polrootsreal(resy,[0,1])
  x0 = 0.35844524659290358614034590968895026895

  polrootsreal(subst(Lx,x,x0),[0,1])
  y0 = 0.24736807727998204016181765766909154018

  my(x=x0,y=y0);(1+x)*(1+x*y)*(1-x*y)^3*(1-y)*(1+y)^3;
  \\ 0.14820813545262897876369097043767859559  <----- value
.........................................................................

/* z = 1, g = 1 */

  L=(1+x)*(1+x*y)*(1-x*y)*(1-x*y*t)^2*(1-y)*(1+y)*(1+y*t)^2

  Tx=factor(deriv(L,x));   \\ it takes some time in PARI use Magma
  Ty=factor(deriv(L,y));
  Tt=factor(deriv(L,t));
  Lx= x^3*y^3*t + 4/5*x^2*y^3*t - 3/5*x^2*y^2 - 2/5*x*y^2 - 3/5*x*y*t - 2/5*y*t + 1/5;
  Ly=Ty[4,1];
  Lt=Tt[9,1];

  T=factor(polresultant(Lx,Lt,t))
  res1t=T[3,1]
  T=factor(polresultant(Ly,Lt,t))
  res2t=T[4,1]

  T=factor(polresultant(res1g,res2g,y))
  resy=T[3,1]

  x0=5/9;

  polrootsreal(subst(res1t,x,x0),[0,1])
  /* no roots */

/* z = 1, g = 1, t = 0
  L=(1+x)*(1-(x*y)^2)*(1-y^2)  <= 2 */

/* z = 1, g = 1, t = 1
  L=(1+x)*(1+x*y)*(1-x*y)^3*(1-y)*(1+y)^3 <= 2 (tested before)  */
.........................................................................

/* g = 1
  (1+x)*(1+x*y)*(1-x*y*z)*(1-x*y*z*t)^2*(1-y)*(1+y*z)*(1+y*z*t)^2 */

  {res1x=
  x - 4/5
  ;}
  
  {res1z=
  x^2*y - 2/5*x*y + 4/5*x - 3/5
  ;}
  
  {res1t=
  x^3*y^2*z + 4/5*x^2*y*z - 4/5*x^2*y - 1/5*x*y*z + 1/5*x*y - 3/5*x + 2/5
  ;}
  
  {Lt=
  x*y*z*t + 1/2*x - 1/2
  ;}

  polrootsreal(res1x,[0,1])
  x0 = 4/5

  polrootsreal(subst(res1z,x,x0),[0,1])
  \\ no roots
.........................................................................
