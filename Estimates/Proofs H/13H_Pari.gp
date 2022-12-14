
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

/* z = 1, t = 1 */
  
  L= (1+x)*(1+x*y)*(1-x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y*a)*(1-y)*(1+y)*(1+2*g*y+y^2)*(1+y*a)
  
  {res1y=  
  x^3 - 61/28*x^2 + 55/28*x - 9/14
  ;}
  
  {res1a=
  x^4*y^4 + 2/7*x^3*y^4 + 5/7*x^3*y^2 - 5/7*x^2*y^2 - 3/7*x + 2/7
  ;}
  
  {res1g=
  x^5*y^5*a + 4/7*x^4*y^5*a - 6/7*x^4*y^4 + 5/7*x^4*y^3*a - 3/7*x^3*y^4 - 3/7*x^3*y^3*a - 4/7*x^3*y^2 
      - 2/7*x^2*y^3*a + 3/7*x^2*y^2 - 3/7*x^2*y*a + 1/7*x*y^2 + 2/7*x - 1/7
  ;}
  
  {Lg=
  x^2*y^2 - x*y^2 - 4*x*y*g - x + 1
  ;}
 
  polrootsreal(res1y,[0,1])
  x0=0.68021174953178485057337263425642396642;
  
  polrootsreal(subst(res1a,x,x0),[0,1])
  y0 = 0.62919611548377409320659100688365835245
  
  polrootsreal(subst(subst(res1g,x,x0),y,y0),[0,1])
  a0 = 0.37359610973533836743840230426411133728
  
  polrootsreal(subst(subst(subst(res1g,x,x0),y,y0),a,a0),[0,1])
  /* no roots */
.........................................................................

/* z = 1, t = 1, a = 0:
  (1+x)*(1+x*y)*(1-x*y)*(1-2*g*x*y+(x*y)^2)*(1-y)*(1+y)*(1+2*g*y+y^2) */

  {res1y=
  x - 5/9
  ;}

  {res1g=
  x^4*y^4 + 1/2*x^3*y^4 + 2/3*x^3*y^2 - 1/2*x^2*y^2 - 1/6*x*y^2 - 1/3*x + 1/6
  ;}

  {Lg=
  x^2*y^2 - x*y^2 - 4*x*y*g - x + 1
  ;}

  polrootsreal(res1y,[0,1])
  x0= 5/9
  
  polrootsreal(subst(res1g,x,x0),[0,1])
  y0 = 0.92338051687663868591100756630627792309
  
  polrootsreal(subst(subst(Lg,x,x0),y,y0),[0,1])
  g0 = 0.11399759467612823282851945263040468186

  my(x=x0,y=y0,g=g0);(1+x)*(1+x*y)*(1-x*y)*(1-2*g*x*y+(x*y)^2)*(1-y)*(1+y)*(1+2*g*y+y^2)
  \\ 0.39944484279605505910808458619888345974  <----- value
.........................................................................

/* z = 1, t = 1, a = 0, g = 0:
  (1+x)*(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-y)*(1+y)*(1+y^2) */

  L=(1+x)*(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-y)*(1+y)*(1+y^2);
  Tx=factor(deriv(L,x));
  Ty=factor(deriv(L,y));
  Lx=Tx[4,1];
  Ly=Ty[3,1];

  T1=factor(polresultant(Lx,Ly,y))
  res1y=T1[3,1]

  polrootsreal(res1y,[0,1])
 /* no roots */
.........................................................................

/* z = 1, t = 1, a = 0, g = 1:
  (1+x)*(1+x*y)*(1-x*y)*(1-2*x*y+(x*y)^2)*(1-y)*(1+y)*(1+2*y+y^2) */

  L=(1+x)*(1+x*y)*(1-x*y)*(1-2*x*y+(x*y)^2)*(1-y)*(1+y)*(1+2*y+y^2);
  Tx=factor(deriv(L,x));
  Ty=factor(deriv(L,y));
  Lx=Tx[4,1];
  Ly=Ty[4,1];

  T1=factor(polresultant(Lx,Ly,y))
  res1y=T1[3,1]

  polrootsreal(res1y,[0,1])
  x0 = 0.35844524659290358614034590968895026895
  
  polrootsreal(subst(Ly,x,x0),[0,1])
  y0 = 0.29951413956028237066885809677779695905
  
  my(x=x0,y=y0); (1+x)*(1+x*y)*(1-x*y)*(1-2*x*y+(x*y)^2)*(1-y)*(1+y)*(1+2*y+y^2)
  \\ 1.6447618467775169468795242068854263610  <----- value
.........................................................................

/* z = 1, t = 1, a = 1:
  (1+x)*(1+x*y)*(1-x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+2*g*y+y^2)*(1+y) */

  {res1y=
  x^7 - 1409/343*x^6 + 20609/2401*x^5 - 26262/2401*x^4 + 21359/2401*x^3 - 1571/343*x^2 + 3315/2401*x -
      450/2401
  ;}

  {res1g=
  x^4*y^4 + 4/7*x^3*y^4 + 1/7*x^3*y^3 + 5/7*x^3*y^2 + 1/7*x^2*y^3 - 2/7*x^2*y^2 + 1/7*x^2*y - 
      1/7*x*y^2 + 1/7*x*y - 2/7*x + 1/7
  ;}

  {Lg=
  x^2*y^2 - x*y^2 - 4*x*y*g - x + 1
  ;}

  polrootsreal(res1y,[0,1])
  x0= 0.61002921924197630535978331353721954447
  
  polrootsreal(subst(res1g,x,x0),[0,1])
  y0 = 0.22361530800201682087887020423477626867
  
  polrootsreal(subst(subst(Lg,x,x0),y,y0),[0,1])
  g0 = 0.69289276648508700595042002289547471665

  my(x=x0,y=y0,g=g0); (1+x)*(1+x*y)*(1-x*y)*(1-2*g*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+2*g*y+y^2)*(1+y)
  \\ 1.7893927819435736927942039186719569639  <----- value
.........................................................................

/* z = 1, t = 1, a = 1, g = 0:
  (1+x)*(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+y^2)*(1+y) */

  L=(1+x)*(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+y^2)*(1+y);
  Tx=factor(deriv(L,x));
  Ty=factor(deriv(L,y));
  Lx=Tx[5,1];
  Ly=Ty[4,1];

  T1=factor(polresultant(Lx,Ly,y))
  res1y=T1[3,1]

  polrootsreal(res1y,[0,1])
  /* no roots */
.........................................................................

/* z = 1, t = 1, a = 1, g = 1:
  (1+x)*(1+x*y)*(1-x*y)*(1-2*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+2*y+y^2)*(1+y) */

  L=(1+x)*(1+x*y)*(1-x*y)*(1-2*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+2*y+y^2)*(1+y);
  Tx=factor(deriv(L,x));
  Ty=factor(deriv(L,y));
  Lx=Tx[4,1];
  Ly=Ty[4,1];

  T1=factor(polresultant(Lx,Ly,y))
  res1y=T1[3,1]

  polrootsreal(res1y,[0,1])
  x0 = 0.61560399916261661720655780480980991771

  polrootsreal(subst(Ly,x,x0),[0,1])
  y0 = 0.17319801367773050963945020285014573842

  my(x=x0,y=y0); (1+x)*(1+x*y)*(1-x*y)*(1-2*x*y+(x*y)^2)*(1-x*y)*(1-y)*(1+y)*(1+2*y+y^2)*(1+y)
  \\ 1.7838690380867688582862509656478683935  <----- value
.........................................................................

/* z = 1, t = 1, g = 0:
  (1+x)*(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-x*y*a)*(1-y)*(1+y)*(1+y^2)*(1+y*a) */  

  {res1y=
  x^4 - 1/2*x^3 + 1/2*x^2 - 1/2*x + 5/6
  ;}
  
  {res1a=
  x^5*y^4 + 1/2*x^4*y^4 - 1/3*x + 1/6
  ;}
  
  {La=
  x*y*a + 1/2*x - 1/2
  ;}
  
  polrootsreal(res1y,[0,1])
  /* no roots */
.........................................................................

/* z = 1, t = 1, g = 1:
  (1+x)*(1+x*y)*(1-x*y)*(1-2*x*y+(x*y)^2)*(1-x*y*a)*(1-y)*(1+y)*(1+2*y+y^2)*(1+y*a) */

  {res1y=
  x^5 - 7/6*x^4 + 85/36*x^3 - 43/24*x^2 + 5/6*x - 25/72
  ;}

  {res1a=
  x^3*y^2 + 1/2*x^2*y^2 + 1/3*x^2*y + 1/3*x*y - 1/3*x + 1/6
  ;}

  {La=
  x*y*a + 1/2*x - 1/2
  ;}

  polrootsreal(res1y,[0,1])
  x0 = 0.65564993996123838841249695254402435326

  polrootsreal(subst(res1a,x,x0),[0,1])
  y0 = 0.12271238894100992184313655816494258613

  polrootsreal(subst(subst(La,x,x0),y,y0),[0,1])
  /* no roots */
.........................................................................

/* z = 1, a = 0, g = 0:
  (1+x)*(1+x*y)*(1-x*y)*(1+(x*y*t)^2)*(1-y)*(1+y)*(1+(y*t)^2) maximized for t=1
  (1+x)*(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-y)*(1+y)*(1+y^2) */

  L= (1+x)*(1+x*y)*(1-x*y)*(1+(x*y)^2)*(1-y)*(1+y)*(1+y^2)
  Tx=factor(deriv(L,x));
  Ty=factor(deriv(L,y));
  Lx=Tx[4,1];
  Ly=Ty[3,1];

  T=factor(polresultant(Lx,Ly,y))
  /* no roots */
.........................................................................

/* z = 1, a = 1:
   (1+x)*(1+x*y)*(1-x*y)*(1-2*g*x*y*t+(x*y*t)^2)*(1-x*y*t)*(1-y)*(1+y)*(1+2*g*y*t+(y*t)^2)*(1+y*t) */

  {res1y=
  x^6 - 197/63*x^5 + 3196/567*x^4 - 3442/567*x^3 + 2551/567*x^2 - 379/189*x + 14/27
  ;}
 
  {res1t=
  x^9*y^6 - 8/63*x^8*y^6 - 53/567*x^7*y^6 - 115/63*x^7*y^4 - 2/567*x^6*y^6 + 806/567*x^6*y^4 + 
      23/189*x^5*y^4 + 625/567*x^5*y^2 - 8/567*x^4*y^4 - 100/63*x^4*y^2 + 35/81*x^3*y^2 - 125/567*x^3 
      + 2/189*x^2*y^2 + 262/567*x^2 - 55/189*x + 4/63
  ;}
  
  {res1g=
  x^5*y^5*t^3 + 4/7*x^4*y^5*t^3 - 6/7*x^4*y^4*t^2 + 5/7*x^4*y^3*t - 3/7*x^3*y^4*t^2 - 5/7*x^3*y^3*t^3 
      + 2/7*x^3*y^3*t - 4/7*x^3*y^2 - 2/7*x^2*y^3*t^3 + 4/7*x^2*y^2*t^2 - 1/7*x^2*y^2 - 3/7*x^2*y*t + 
      1/7*x*y^2*t^2 + 2/7*x - 1/7
  ;}

  {Lg=
  x^2*y^2*t^2 - x*y^2*t^2 - 4*x*y*t*g - x + 1
  ;}

  polrootsreal(res1y,[0,1])
  /* no roots */
.........................................................................

/* z = 1, a = 1, g = 0:
  (1+x)*(1+x*y)*(1-x*y)*(1+(x*y*t)^2)*(1-x*y*t)*(1-y)*(1+y)*(1+(y*t)^2)*(1+y*t) */

  {res1y=
  x^14 - 29/18*x^13 + 116/27*x^12 - 535/108*x^11 + 899/108*x^10 - 439/54*x^9 + 587/54*x^8 - 710/81*x^7
      + 749/81*x^6 - 979/162*x^5 + 1297/243*x^4 - 2677/972*x^3 + 1865/972*x^2 - 275/486*x + 125/486
  ;}

  {res1t=
  x^19*y^10 + 25/18*x^18*y^10 + 43/54*x^17*y^10 - 8/3*x^17*y^8 + 13/108*x^16*y^10 - 155/54*x^16*y^8 + 
      1/54*x^15*y^10 - 40/27*x^15*y^8 + 25/9*x^15*y^6 + 53/324*x^14*y^10 - 11/108*x^14*y^8 + 
      323/162*x^14*y^6 + 53/486*x^13*y^10 - 118/243*x^13*y^8 + 247/243*x^13*y^6 - 38/27*x^13*y^4 + 
      1/36*x^12*y^10 - 749/972*x^12*y^8 - 22/243*x^12*y^6 - 205/486*x^12*y^4 - 1/486*x^11*y^10 - 
      22/81*x^11*y^8 + 23/27*x^11*y^6 - 89/243*x^11*y^4 + 28/81*x^11*y^2 + 1/972*x^10*y^10 - 
      17/972*x^10*y^8 + 161/486*x^10*y^6 + 11/81*x^10*y^4 - 16/243*x^10*y^2 + 4/243*x^9*y^8 - 
      25/81*x^9*y^6 - 140/243*x^9*y^4 + 17/162*x^9*y^2 - 8/243*x^9 - 1/486*x^8*y^10 - 5/324*x^8*y^8 - 
      91/486*x^8*y^6 + 79/162*x^8*y^4 - 53/972*x^8*y^2 + 2/81*x^8 + 5/243*x^7*y^8 + 130/243*x^7*y^4 + 
      89/486*x^7*y^2 - 5/243*x^7 + 1/54*x^6*y^8 - 2/243*x^6*y^6 + 35/162*x^6*y^4 - 41/108*x^6*y^2 + 
      7/972*x^6 - 32/243*x^5*y^6 + 2/243*x^5*y^4 - 79/486*x^5*y^2 - 2/81*x^5 - 1/18*x^4*y^6 + 
      43/243*x^4*y^4 - 85/972*x^4*y^2 + 67/972*x^4 + 19/81*x^3*y^4 - 1/162*x^3*y^2 + 1/18*x^2*y^4 - 
      43/324*x^2*y^2 + 5/324*x^2 - 2/27*x*y^2 + 1/36
  ;}

  {Lt=
  x^3*y^5*t^5 + 5/6*x^3*y^4*t^4 + 2/3*x^3*y^3*t^3 + 1/2*x^3*y^2*t^2 - 5/6*x^2*y^4*t^4 - 
      2/3*x^2*y^3*t^3 - 1/2*x^2*y^2*t^2 - 1/3*x^2*y*t + 2/3*x*y^3*t^3 + 1/2*x*y^2*t^2 + 1/3*x*y*t + 
      1/6*x - 1/2*y^2*t^2 - 1/3*y*t - 1/6
  ;}

  polrootsreal(res1y,[0,1])
  /* no roots */
.........................................................................

/* t = 0:
  (1+x)*(1+x*y)*(1-x*y*z)*(1-y)*(1+y*z) */

  L=(1+x)*(1+x*y)*(1-x*y*z)*(1-y)*(1+y*z);
  Tx=factor(deriv(L,x));
  Ty=factor(deriv(L,y));
  Tz=factor(deriv(L,z));
  Lx=Tx[3,1];
  Ly=Ty[2,1];
  Lz=Tz[5,1];

  T1=factor(polresultant(Lx,Lz,z))
  T2=factor(polresultant(Ly,Lz,z))
  res1z=T1[3,1]
  res2z=T2[4,1]

  T11=factor(polresultant(res1z,res2z,y))
  res1y=T11[3,1]

  polrootsreal(res1y,[0,1])
  x0 = 2/3

  polrootsreal(subst(res1z,x,x0),[0,1])
  /* no roots */
.........................................................................

/* t = 1:
  (1+x)*(1+x*y)*(1-x*y*z)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z*a)*(1-y)*(1+y*z)*(1+2*g*y*z+(y*z)^2)*(1+y*z*a) */

  {res1y=
  x^3 - 52/21*x^2 + 1123/567*x - 14/27
  ;}

  {res1z=
  x^6*y^3 - 82/63*x^5*y^3 + 152/63*x^5*y^2 + 31/63*x^4*y^3 - 751/189*x^4*y^2 + 157/81*x^4*y - 
      4/63*x^3*y^3 + 358/189*x^3*y^2 - 2192/567*x^3*y + 14/27*x^3 - 55/189*x^2*y^2 + 1319/567*x^2*y - 
      689/567*x^2 - 250/567*x*y + 512/567*x - 125/567
  ;}

  {res1a=
  x^4*y^4*z^3 + 2/7*x^3*y^4*z^3 + 6/7*x^3*y^3*z^3 - 6/7*x^3*y^3*z^2 + 5/7*x^3*y^2*z + 1/7*x^2*y^3*z^3 
      - 1/7*x^2*y^3*z^2 - 5/7*x^2*y^2*z^2 + 4/7*x^2*y*z - 4/7*x^2*y - 1/7*x*y*z + 1/7*x*y - 3/7*x + 
      2/7
  ;}

  {res1g=
  x^5*y^5*z^4*a + 4/7*x^4*y^5*z^4*a + 6/7*x^4*y^4*z^4*a - 6/7*x^4*y^4*z^3*a - 6/7*x^4*y^4*z^3 + 
      5/7*x^4*y^3*z^2*a + 3/7*x^3*y^4*z^4*a - 3/7*x^3*y^4*z^3*a - 3/7*x^3*y^4*z^3 - 5/7*x^3*y^3*z^3*a 
      - 5/7*x^3*y^3*z^3 + 2/7*x^3*y^3*z^2*a + 5/7*x^3*y^3*z^2 + 4/7*x^3*y^2*z^2*a - 4/7*x^3*y^2*z*a - 
      4/7*x^3*y^2*z - 2/7*x^2*y^3*z^3*a - 2/7*x^2*y^3*z^3 + 2/7*x^2*y^3*z^2 + 1/7*x^2*y^2*z^2*a + 
      4/7*x^2*y^2*z^2 - 1/7*x^2*y^2*z*a - 1/7*x^2*y^2*z - 3/7*x^2*y*z*a - 3/7*x^2*y*z + 3/7*x^2*y + 
      1/7*x*y^2*z^2 + 2/7*x - 1/7
  ;}

  {Lg=
  x^2*y^2*z^2 - x*y^2*z^2 - 4*x*y*z*g - x + 1
  ;}

  polrootsreal(res1y,[0,1])
  /* no roots */
.........................................................................

/* t = 1, a = 0:
  (1+x)*(1+x*y)*(1-x*y*z)*(1-2*g*x*y*z+(x*y*z)^2)*(1-y)*(1+y*z)*(1+2*g*y*z+(y*z)^2) */

  {res1y=
  x^3 - 689/294*x^2 + 256/147*x - 125/294
  ;}

  {res1z=
  x^6*y^3 - 85/98*x^5*y^3 + 95/42*x^5*y^2 + 6/49*x^4*y^3 - 424/147*x^4*y^2 + 250/147*x^4*y - 
      1/98*x^3*y^3 + 229/294*x^3*y^2 - 281/98*x^3*y + 125/294*x^3 - 11/147*x^2*y^2 + 185/147*x^2*y - 
      131/147*x^2 - 17/98*x*y + 55/98*x - 6/49
  ;}

  {res1g=
  x^4*y^4*z^3 + 1/2*x^3*y^4*z^3 + 5/6*x^3*y^3*z^3 - 5/6*x^3*y^3*z^2 + 2/3*x^3*y^2*z + 1/3*x^2*y^3*z^3 
      - 1/3*x^2*y^3*z^2 - 2/3*x^2*y^2*z^2 + 1/6*x^2*y^2*z + 1/2*x^2*y*z - 1/2*x^2*y - 1/6*x*y^2*z^2 - 
      1/3*x + 1/6
  ;}

  {Lg=
  x^2*y^2*z^2 - x*y^2*z^2 - 4*x*y*z*g - x + 1
  ;}

  polrootsreal(res1y,[0,1])
  /* no roots */
.........................................................................

/* t = 1, a = 0, g = 0:
  (1+x)*(1+x*y)*(1-x*y*z)*(1+(x*y*z)^2)*(1-y)*(1+y*z)*(1+(y*z)^2) */

  {res1y=
  x^11 - 12/5*x^10 + 33/10*x^9 - 506/125*x^8 + 257/50*x^7 - 744/125*x^6 + 49/10*x^5 - 458/125*x^4 + 
      683/250*x^3 - 242/125*x^2 + 4/5*x - 16/125
  ;}

  {res1z=
  x^14*y^5 + 4/5*x^13*y^5 + 17/5*x^13*y^4 + 23/50*x^12*y^5 + 47/25*x^12*y^4 + 457/100*x^12*y^3 - 
      1/125*x^11*y^5 + 156/125*x^11*y^4 + 159/125*x^11*y^3 + 379/125*x^11*y^2 + 13/50*x^10*y^5 - 
      9/50*x^10*y^4 + 717/500*x^10*y^3 - 43/500*x^10*y^2 + 124/125*x^10*y + 1/25*x^9*y^5 + 
      29/25*x^9*y^4 - 109/250*x^9*y^3 + 241/250*x^9*y^2 - 46/125*x^9*y + 16/125*x^9 - 19/250*x^8*y^5 -
      23/50*x^8*y^4 + 477/250*x^8*y^3 - 207/500*x^8*y^2 + 101/250*x^8*y - 12/125*x^8 - 7/125*x^7*y^5 -
      79/125*x^7*y^4 - 249/125*x^7*y^3 + 189/125*x^7*y^2 - 22/125*x^7*y + 2/25*x^7 - 1/250*x^6*y^5 - 
      87/250*x^6*y^4 - 351/250*x^6*y^3 - 131/50*x^6*y^2 + 149/250*x^6*y - 7/250*x^6 - 1/25*x^5*y^5 - 
      2/125*x^5*y^4 - 92/125*x^5*y^3 - 144/125*x^5*y^2 - 176/125*x^5*y + 12/125*x^5 - 7/125*x^4*y^5 - 
      73/250*x^4*y^4 - 3/100*x^4*y^3 - 7/10*x^4*y^2 - 79/250*x^4*y - 67/250*x^4 - 2/125*x^3*y^5 - 
      39/125*x^3*y^4 - 98/125*x^3*y^3 - 4/125*x^3*y^2 - 8/25*x^3*y - 9/125*x^2*y^4 - 311/500*x^2*y^3 -
      479/500*x^2*y^2 - 3/250*x^2*y - 3/50*x^2 - 27/250*x*y^3 - 129/250*x*y^2 - 66/125*x*y - 
      27/500*y^2 - 18/125*y - 27/250
  ;}

  {Lz=
  x^3*y^5*z^5 + 5/6*x^3*y^4*z^4 + 2/3*x^3*y^3*z^3 + 1/2*x^3*y^2*z^2 - 5/6*x^2*y^4*z^4 - 
      2/3*x^2*y^3*z^3 - 1/2*x^2*y^2*z^2 - 1/3*x^2*y*z + 2/3*x*y^3*z^3 + 1/2*x*y^2*z^2 + 1/3*x*y*z + 
      1/6*x - 1/2*y^2*z^2 - 1/3*y*z - 1/6
  ;}

  polrootsreal(res1y,[0,1])
  /* no roots */
.........................................................................

/* t = 1, a = 0, g = 1:
  (1+x)*(1+x*y)*(1-x*y*z)*(1-2*x*y*z+(x*y*z)^2)*(1-y)*(1+y*z)*(1+2*y*z+(y*z)^2) */

  {res1y=
  x - 4/5
  ;}

  {res1z=
  x^2*y - 2/5*x*y + 4/5*x - 3/5
  ;}

  {Lz=
  x*y*z + 1/2*x - 1/2
  ;}

  polrootsreal(res1y,[0,1])
  x0 = 4/5;
  
  polrootsreal(subst(res1z,x,x0),[0,1])
  /* no roots */
.........................................................................

/* t = 1, a = 1:
  (1+x)*(1+x*y)*(1-x*y*z)*(1-2*g*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+y*z)*(1+2*g*y*z+(y*z)^2)*(1+y*z) */

  {res1y=
  x - 27/28
  ;}

  {res1z=
  x^6*y^3 - 13/14*x^5*y^3 + 33/14*x^5*y^2 + 23/112*x^4*y^3 - 165/56*x^4*y^2 + 207/112*x^4*y - 
      1/28*x^3*y^3 + 57/56*x^3*y^2 - 81/28*x^3*y + 27/56*x^3 - 3/16*x^2*y^2 + 81/56*x^2*y - 
      101/112*x^2 - 9/28*x*y + 17/28*x - 5/28
  ;}

  {res1g=
  x^4*y^4*z^3 + 4/7*x^3*y^4*z^3 + 6/7*x^3*y^3*z^3 - 5/7*x^3*y^3*z^2 + 5/7*x^3*y^2*z + 3/7*x^2*y^3*z^3 
      - 2/7*x^2*y^3*z^2 - 4/7*x^2*y^2*z^2 + 2/7*x^2*y^2*z + 4/7*x^2*y*z - 3/7*x^2*y - 1/7*x*y^2*z^2 + 
      1/7*x*y*z - 2/7*x + 1/7
  ;}

  {Lg=
  x^2*y^2*z^2 - x*y^2*z^2 - 4*x*y*z*g - x + 1
  ;}

  polrootsreal(res1y,[0,1])
  x0 = 27/28;
  
  polrootsreal(subst(res1z,x,x0),[0,1])
  /* no roots */
.........................................................................

/* t = 1, a = 1, g = 0:
  (1+x)*(1+x*y)*(1-x*y*z)*(1+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+y*z)*(1+(y*z)^2)*(1+y*z) */

  {res1y=
  x^11 - 5/2*x^10 + 3*x^9 - 115/54*x^8 + 1481/432*x^7 - 4849/864*x^6 + 919/216*x^5 - 1639/864*x^4 + 
      869/432*x^3 - 1955/864*x^2 + 25/27*x - 125/864
  ;}

  {res1z=
  x^14*y^5 + 1/2*x^13*y^5 + 7/2*x^13*y^4 + 5/6*x^12*y^4 + 29/6*x^12*y^3 + 8/27*x^11*y^5 + 
      1/36*x^11*y^4 - 5/18*x^11*y^3 + 355/108*x^11*y^2 + 71/144*x^10*y^5 + 439/216*x^10*y^4 + 
      13/27*x^10*y^3 - 31/24*x^10*y^2 + 475/432*x^10*y - 391/864*x^9*y^5 + 1733/864*x^9*y^4 + 
      2003/432*x^9*y^3 + 149/144*x^9*y^2 - 245/288*x^9*y + 125/864*x^9 - 4/27*x^8*y^5 - 
      445/144*x^8*y^4 + 101/36*x^8*y^3 + 1033/216*x^8*y^2 + 85/108*x^8*y - 25/144*x^8 + 13/32*x^7*y^5 
      - 569/864*x^7*y^4 - 1057/144*x^7*y^3 + 661/432*x^7*y^2 + 1991/864*x^7*y + 175/864*x^7 + 
      85/432*x^6*y^5 + 467/216*x^6*y^4 - 149/216*x^6*y^3 - 95/12*x^6*y^2 + 23/144*x^6*y + 91/216*x^6 -
      13/96*x^5*y^5 + 79/96*x^5*y^4 + 1961/432*x^5*y^3 + 133/432*x^5*y^2 - 3437/864*x^5*y - 65/864*x^5
      - 1/8*x^4*y^5 - 43/48*x^4*y^4 + 7/6*x^4*y^3 + 1001/216*x^4*y^2 + 163/216*x^4*y - 109/144*x^4 - 
      1/32*x^3*y^5 - 59/96*x^3*y^4 - 961/432*x^3*y^3 + 235/432*x^3*y^2 + 655/288*x^3*y + 79/288*x^3 - 
      1/8*x^2*y^4 - 79/72*x^2*y^3 - 187/72*x^2*y^2 - 19/216*x^2*y + 23/54*x^2 - 1/6*x*y^3 - 5/6*x*y^2 
      - 77/54*x*y - 5/54*x - 2/27*y^2 - 2/9*y - 8/27
  ;}

  {Lz=
  x^3*y^5*z^5 + 3/4*x^3*y^4*z^4 + 3/4*x^3*y^3*z^3 + 1/2*x^3*y^2*z^2 - 3/4*x^2*y^4*z^4 - 
      1/2*x^2*y^3*z^3 - 1/2*x^2*y^2*z^2 - 1/4*x^2*y*z + 3/4*x*y^3*z^3 + 1/2*x*y^2*z^2 + 1/2*x*y*z + 
      1/4*x - 1/2*y^2*z^2 - 1/4*y*z - 1/4
  ;}

  polrootsreal(res1y,[0,1])
  x0 = 0.96193505139403980459631128904862578003;
  
  polrootsreal(subst(res1z,x,x0),[0,1])
  /* no roots */
.........................................................................

/* t = 1, a = 1, g = 1:
  (1+x)*(1+x*y)*(1-x*y*z)*(1-2*x*y*z+(x*y*z)^2)*(1-x*y*z)*(1-y)*(1+y*z)*(1+2*y*z+(y*z)^2)*(1+y*z) */

  {res1y=
  x^11 - 5/2*x^10 + 3*x^9 - 115/54*x^8 + 1481/432*x^7 - 4849/864*x^6 + 919/216*x^5 - 1639/864*x^4 + 
      869/432*x^3 - 1955/864*x^2 + 25/27*x - 125/864
  ;}

  {res1z=
  x^14*y^5 + 1/2*x^13*y^5 + 7/2*x^13*y^4 + 5/6*x^12*y^4 + 29/6*x^12*y^3 + 8/27*x^11*y^5 + 
      1/36*x^11*y^4 - 5/18*x^11*y^3 + 355/108*x^11*y^2 + 71/144*x^10*y^5 + 439/216*x^10*y^4 + 
      13/27*x^10*y^3 - 31/24*x^10*y^2 + 475/432*x^10*y - 391/864*x^9*y^5 + 1733/864*x^9*y^4 + 
      2003/432*x^9*y^3 + 149/144*x^9*y^2 - 245/288*x^9*y + 125/864*x^9 - 4/27*x^8*y^5 - 
      445/144*x^8*y^4 + 101/36*x^8*y^3 + 1033/216*x^8*y^2 + 85/108*x^8*y - 25/144*x^8 + 13/32*x^7*y^5 
      - 569/864*x^7*y^4 - 1057/144*x^7*y^3 + 661/432*x^7*y^2 + 1991/864*x^7*y + 175/864*x^7 + 
      85/432*x^6*y^5 + 467/216*x^6*y^4 - 149/216*x^6*y^3 - 95/12*x^6*y^2 + 23/144*x^6*y + 91/216*x^6 -
      13/96*x^5*y^5 + 79/96*x^5*y^4 + 1961/432*x^5*y^3 + 133/432*x^5*y^2 - 3437/864*x^5*y - 65/864*x^5
      - 1/8*x^4*y^5 - 43/48*x^4*y^4 + 7/6*x^4*y^3 + 1001/216*x^4*y^2 + 163/216*x^4*y - 109/144*x^4 - 
      1/32*x^3*y^5 - 59/96*x^3*y^4 - 961/432*x^3*y^3 + 235/432*x^3*y^2 + 655/288*x^3*y + 79/288*x^3 - 
      1/8*x^2*y^4 - 79/72*x^2*y^3 - 187/72*x^2*y^2 - 19/216*x^2*y + 23/54*x^2 - 1/6*x*y^3 - 5/6*x*y^2 
      - 77/54*x*y - 5/54*x - 2/27*y^2 - 2/9*y - 8/27
  ;}

  {Lz=
  x^3*y^5*z^5 + 3/4*x^3*y^4*z^4 + 3/4*x^3*y^3*z^3 + 1/2*x^3*y^2*z^2 - 3/4*x^2*y^4*z^4 - 
      1/2*x^2*y^3*z^3 - 1/2*x^2*y^2*z^2 - 1/4*x^2*y*z + 3/4*x*y^3*z^3 + 1/2*x*y^2*z^2 + 1/2*x*y*z + 
      1/4*x - 1/2*y^2*z^2 - 1/4*y*z - 1/4
  ;}

  polrootsreal(res1y,[0,1])
  x0 = 0.96193505139403980459631128904862578003;
  
  polrootsreal(subst(res1z,x,x0),[0,1])
  /* no roots */
.........................................................................

/* a = 0, g = 1:
  (1+x)*(1+x*y)*(1-x*y*z)*(1-2*x*y*z*t+(x*y*z*t)^2)*(1-y)*(1+y*z)*(1+2*y*z*t+(y*z*t)^2) */

  {res1y=
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

  polrootsreal(res1y,[0,1])
  x0 = 4/5
  
  polrootsreal(subst(res1z,x,x0),[0,1])
  /* no roots */

.........................................................................
/* a = 1:
  (1+x)*(1+x*y)*(1-x*y*z)*(1-2*g*x*y*z*t+(x*y*z*t)^2)*(1-x*y*z*t)*(1-y)*(1+y*z)*(1+2*g*y*z*t+(y*z*t)^2)*(1+y*z*t) */

  {res1y=
  x^3 - 52/21*x^2 + 1123/567*x - 14/27
  ;}

  {res1z=
  x^6*y^3 - 82/63*x^5*y^3 + 152/63*x^5*y^2 + 31/63*x^4*y^3 - 751/189*x^4*y^2 + 157/81*x^4*y - 
      4/63*x^3*y^3 + 358/189*x^3*y^2 - 2192/567*x^3*y + 14/27*x^3 - 55/189*x^2*y^2 + 1319/567*x^2*y - 
      689/567*x^2 - 250/567*x*y + 512/567*x - 125/567
  ;}

  {res1t=
  x^9*y^6*z^3 - 8/63*x^8*y^6*z^3 + 152/63*x^8*y^5*z^3 - 152/63*x^8*y^5*z^2 - 53/567*x^7*y^6*z^3 - 
      73/81*x^7*y^5*z^3 + 73/81*x^7*y^5*z^2 + 157/81*x^7*y^4*z^3 - 3233/567*x^7*y^4*z^2 + 
      157/81*x^7*y^4*z - 2/567*x^6*y^6*z^3 - 38/189*x^6*y^5*z^3 + 38/189*x^6*y^5*z^2 - 
      682/567*x^6*y^4*z^3 + 310/81*x^6*y^4*z^2 - 682/567*x^6*y^4*z + 14/27*x^6*y^3*z^3 - 
      2542/567*x^6*y^3*z^2 + 2542/567*x^6*y^3*z - 14/27*x^6*y^3 + 1/567*x^5*y^5*z^3 - 
      1/567*x^5*y^5*z^2 - 13/567*x^5*y^4*z^3 + 95/567*x^5*y^4*z^2 - 13/567*x^5*y^4*z - 
      85/189*x^5*y^3*z^3 + 353/81*x^5*y^3*z^2 - 353/81*x^5*y^3*z + 85/189*x^5*y^3 - 95/81*x^5*y^2*z^2 
      + 1955/567*x^5*y^2*z - 95/81*x^5*y^2 + 4/567*x^4*y^4*z^3 - 16/567*x^4*y^4*z^2 + 4/567*x^4*y^4*z 
      + 4/63*x^4*y^3*z^3 - 284/567*x^4*y^3*z^2 + 284/567*x^4*y^3*z - 4/63*x^4*y^3 + 
      848/567*x^4*y^2*z^2 - 2596/567*x^4*y^2*z + 848/567*x^4*y^2 + 500/567*x^4*y*z - 500/567*x^4*y - 
      1/189*x^3*y^3*z^3 - 5/567*x^3*y^3*z^2 + 5/567*x^3*y^3*z + 1/189*x^3*y^3 - 229/567*x^3*y^2*z^2 + 
      703/567*x^3*y^2*z - 229/567*x^3*y^2 - 281/189*x^3*y*z + 281/189*x^3*y - 125/567*x^3 + 
      22/567*x^2*y^2*z^2 - 38/567*x^2*y^2*z + 22/567*x^2*y^2 + 370/567*x^2*y*z - 370/567*x^2*y + 
      262/567*x^2 - 17/189*x*y*z + 17/189*x*y - 55/189*x + 4/63
  ;}

  {res1g=
  x^5*y^5*z^4*t^3 + 4/7*x^4*y^5*z^4*t^3 + 6/7*x^4*y^4*z^4*t^3 - 6/7*x^4*y^4*z^3*t^3 - 
      6/7*x^4*y^4*z^3*t^2 + 5/7*x^4*y^3*z^2*t + 3/7*x^3*y^4*z^4*t^3 - 3/7*x^3*y^4*z^3*t^3 - 
      3/7*x^3*y^4*z^3*t^2 - 5/7*x^3*y^3*z^3*t^3 - 5/7*x^3*y^3*z^3*t^2 + 5/7*x^3*y^3*z^2*t^2 + 
      2/7*x^3*y^3*z^2*t + 4/7*x^3*y^2*z^2*t - 4/7*x^3*y^2*z*t - 4/7*x^3*y^2*z - 2/7*x^2*y^3*z^3*t^3 - 
      2/7*x^2*y^3*z^3*t^2 + 2/7*x^2*y^3*z^2*t^2 + 4/7*x^2*y^2*z^2*t^2 + 1/7*x^2*y^2*z^2*t - 
      1/7*x^2*y^2*z*t - 1/7*x^2*y^2*z - 3/7*x^2*y*z*t - 3/7*x^2*y*z + 3/7*x^2*y + 1/7*x*y^2*z^2*t^2 + 
      2/7*x - 1/7
  ;}

  {Lg=
  x^2*y^2*z^2*t^2 - x*y^2*z^2*t^2 - 4*x*y*z*t*g - x + 1
  ;}

  polrootsreal(res1y,[0,1])
  /* no roots */
.........................................................................

/* a = 1, g = 0:
  (1+x)*(1+x*y)*(1-x*y*z)*(1+(x*y*z*t)^2)*(1-x*y*z*t)*(1-y)*(1+y*z)*(1+(y*z*t)^2)*(1+y*z*t) */

  {res1y=
  x^11 - 59/18*x^10 + 311/54*x^9 - 2621/324*x^8 + 5077/486*x^7 - 11485/972*x^6 + 1735/162*x^5 - 
      299/36*x^4 + 2867/486*x^3 - 3475/972*x^2 + 350/243*x - 125/486
  ;}

  {res1z=
  x^16*y^5 - 17/18*x^15*y^5 + 23/6*x^15*y^4 + 7/6*x^14*y^5 - 251/54*x^14*y^4 + 211/36*x^14*y^3 - 
      449/324*x^13*y^5 + 1777/324*x^13*y^4 - 5665/648*x^13*y^3 + 965/216*x^13*y^2 + 29/18*x^12*y^5 - 
      3199/486*x^12*y^4 + 5033/486*x^12*y^3 - 7739/972*x^12*y^2 + 275/162*x^12*y - 5/4*x^11*y^5 + 
      7693/972*x^11*y^4 - 24239/1944*x^11*y^3 + 18965/1944*x^11*y^2 - 860/243*x^11*y + 125/486*x^11 + 
      127/162*x^10*y^5 - 3325/486*x^10*y^4 + 1241/81*x^10*y^3 - 106/9*x^10*y^2 + 370/81*x^10*y - 
      50/81*x^10 - 193/324*x^9*y^5 + 4237/972*x^9*y^4 - 175/12*x^9*y^3 + 14269/972*x^9*y^2 - 
      599/108*x^9*y + 275/324*x^9 + 59/162*x^8*y^5 - 1609/486*x^8*y^4 + 4711/486*x^8*y^3 - 
      3680/243*x^8*y^2 + 1697/243*x^8*y - 253/243*x^8 - 17/108*x^7*y^5 + 2095/972*x^7*y^4 - 
      7087/972*x^7*y^3 + 10417/972*x^7*y^2 - 7471/972*x^7*y + 1285/972*x^7 - 29/27*x^6*y^4 + 
      539/108*x^6*y^3 - 3877/486*x^6*y^2 + 1423/243*x^6*y - 124/81*x^6 + 5/54*x^5*y^4 - 
      5447/1944*x^5*y^3 + 1225/216*x^5*y^2 - 4223/972*x^5*y + 1225/972*x^5 + 38/81*x^4*y^3 - 
      3413/972*x^4*y^2 + 1543/486*x^4*y - 229/243*x^4 - 1/72*x^3*y^3 + 1697/1944*x^3*y^2 - 
      2063/972*x^3*y + 683/972*x^3 - 1/18*x^2*y^2 + 19/27*x^2*y - 121/243*x^2 - 2/27*x*y + 50/243*x - 
      8/243
  ;}

  {res1t=
  x^19*y^10*z^5 + 25/18*x^18*y^10*z^5 + 23/6*x^18*y^9*z^5 - 23/6*x^18*y^9*z^4 + 43/54*x^17*y^10*z^5 + 
      265/54*x^17*y^9*z^5 - 265/54*x^17*y^9*z^4 + 211/36*x^17*y^8*z^5 - 259/18*x^17*y^8*z^4 + 
      211/36*x^17*y^8*z^3 + 13/108*x^16*y^10*z^5 + 295/108*x^16*y^9*z^5 - 295/108*x^16*y^9*z^4 + 
      4433/648*x^16*y^8*z^5 - 5363/324*x^16*y^8*z^4 + 4433/648*x^16*y^8*z^3 + 965/216*x^16*y^7*z^5 - 
      1549/72*x^16*y^7*z^4 + 1549/72*x^16*y^7*z^3 - 965/216*x^16*y^7*z^2 + 1/54*x^15*y^10*z^5 + 
      19/54*x^15*y^9*z^5 - 19/54*x^15*y^9*z^4 + 1811/486*x^15*y^8*z^5 - 2171/243*x^15*y^8*z^4 + 
      1811/486*x^15*y^8*z^3 + 2279/486*x^15*y^7*z^5 - 1766/81*x^15*y^7*z^4 + 1766/81*x^15*y^7*z^3 - 
      2279/486*x^15*y^7*z^2 + 275/162*x^15*y^6*z^5 - 1297/81*x^15*y^6*z^4 + 848/27*x^15*y^6*z^3 - 
      1297/81*x^15*y^6*z^2 + 275/162*x^15*y^6*z + 53/324*x^14*y^10*z^5 + 281/972*x^14*y^9*z^5 - 
      281/972*x^14*y^9*z^4 + 739/1944*x^14*y^8*z^5 - 419/486*x^14*y^8*z^4 + 739/1944*x^14*y^8*z^3 + 
      4925/1944*x^14*y^7*z^5 - 22511/1944*x^14*y^7*z^4 + 22511/1944*x^14*y^7*z^3 - 
      4925/1944*x^14*y^7*z^2 + 85/54*x^14*y^6*z^5 - 9013/648*x^14*y^6*z^4 + 8639/324*x^14*y^6*z^3 - 
      9013/648*x^14*y^6*z^2 + 85/54*x^14*y^6*z + 125/486*x^14*y^5*z^5 - 5765/972*x^14*y^5*z^4 + 
      22147/972*x^14*y^5*z^3 - 22147/972*x^14*y^5*z^2 + 5765/972*x^14*y^5*z - 125/486*x^14*y^5 + 
      53/486*x^13*y^10*z^5 + 193/243*x^13*y^9*z^5 - 193/243*x^13*y^9*z^4 + 181/243*x^13*y^8*z^5 - 
      160/81*x^13*y^8*z^4 + 181/243*x^13*y^8*z^3 + 19/108*x^13*y^7*z^5 - 623/972*x^13*y^7*z^4 + 
      623/972*x^13*y^7*z^3 - 19/108*x^13*y^7*z^2 + 419/486*x^13*y^6*z^5 - 7261/972*x^13*y^6*z^4 + 
      6917/486*x^13*y^6*z^3 - 7261/972*x^13*y^6*z^2 + 419/486*x^13*y^6*z + 50/243*x^13*y^5*z^5 - 
      1031/243*x^13*y^5*z^4 + 7447/486*x^13*y^5*z^3 - 7447/486*x^13*y^5*z^2 + 1031/243*x^13*y^5*z - 
      50/243*x^13*y^5 - 425/486*x^13*y^4*z^4 + 7981/972*x^13*y^4*z^3 - 2605/162*x^13*y^4*z^2 + 
      7981/972*x^13*y^4*z - 425/486*x^13*y^4 + 1/36*x^12*y^10*z^5 + 397/972*x^12*y^9*z^5 - 
      397/972*x^12*y^9*z^4 + 409/324*x^12*y^8*z^5 - 3203/972*x^12*y^8*z^4 + 409/324*x^12*y^8*z^3 + 
      757/972*x^12*y^7*z^5 - 4043/972*x^12*y^7*z^4 + 4043/972*x^12*y^7*z^3 - 757/972*x^12*y^7*z^2 + 
      25/972*x^12*y^6*z^5 - 35/1944*x^12*y^6*z^4 - 103/972*x^12*y^6*z^3 - 35/1944*x^12*y^6*z^2 + 
      25/972*x^12*y^6*z + 115/972*x^12*y^5*z^5 - 1177/486*x^12*y^5*z^4 + 8503/972*x^12*y^5*z^3 - 
      8503/972*x^12*y^5*z^2 + 1177/486*x^12*y^5*z - 115/972*x^12*y^5 - 235/486*x^12*y^4*z^4 + 
      7705/1944*x^12*y^4*z^3 - 7175/972*x^12*y^4*z^2 + 7705/1944*x^12*y^4*z - 235/486*x^12*y^4 + 
      2285/1944*x^12*y^3*z^3 - 45/8*x^12*y^3*z^2 + 45/8*x^12*y^3*z - 2285/1944*x^12*y^3 - 
      1/486*x^11*y^10*z^5 + 19/243*x^11*y^9*z^5 - 19/243*x^11*y^9*z^4 + 227/486*x^11*y^8*z^5 - 
      293/243*x^11*y^8*z^4 + 227/486*x^11*y^8*z^3 + 137/162*x^11*y^7*z^5 - 2143/486*x^11*y^7*z^4 + 
      2143/486*x^11*y^7*z^3 - 137/162*x^11*y^7*z^2 + 10/27*x^11*y^6*z^5 - 635/162*x^11*y^6*z^4 + 
      644/81*x^11*y^6*z^3 - 635/162*x^11*y^6*z^2 + 10/27*x^11*y^6*z - 1/486*x^11*y^5*z^5 + 
      38/243*x^11*y^5*z^4 - 194/243*x^11*y^5*z^3 + 194/243*x^11*y^5*z^2 - 38/243*x^11*y^5*z + 
      1/486*x^11*y^5 - 26/81*x^11*y^4*z^4 + 223/81*x^11*y^4*z^3 - 1271/243*x^11*y^4*z^2 + 
      223/81*x^11*y^4*z - 26/81*x^11*y^4 + 53/162*x^11*y^3*z^3 - 103/81*x^11*y^3*z^2 + 
      103/81*x^11*y^3*z - 53/162*x^11*y^3 - 379/486*x^11*y^2*z^2 + 463/243*x^11*y^2*z - 
      379/486*x^11*y^2 + 1/972*x^10*y^10*z^5 - 13/972*x^10*y^9*z^5 + 13/972*x^10*y^9*z^4 + 
      13/324*x^10*y^8*z^5 - 95/972*x^10*y^8*z^4 + 13/324*x^10*y^8*z^3 + 17/108*x^10*y^7*z^5 - 
      229/324*x^10*y^7*z^4 + 229/324*x^10*y^7*z^3 - 17/108*x^10*y^7*z^2 + 215/972*x^10*y^6*z^5 - 
      2111/972*x^10*y^6*z^4 + 2057/486*x^10*y^6*z^3 - 2111/972*x^10*y^6*z^2 + 215/972*x^10*y^6*z + 
      65/972*x^10*y^5*z^5 - 1691/972*x^10*y^5*z^4 + 3415/486*x^10*y^5*z^3 - 3415/486*x^10*y^5*z^2 + 
      1691/972*x^10*y^5*z - 65/972*x^10*y^5 + 5/108*x^10*y^4*z^4 - 349/648*x^10*y^4*z^3 + 
      121/108*x^10*y^4*z^2 - 349/648*x^10*y^4*z + 5/108*x^10*y^4 + 239/648*x^10*y^3*z^3 - 
      3271/1944*x^10*y^3*z^2 + 3271/1944*x^10*y^3*z - 239/648*x^10*y^3 + 43/1944*x^10*y^2*z^2 - 
      107/972*x^10*y^2*z + 43/1944*x^10*y^2 + 62/243*x^10*y*z - 62/243*x^10*y + 5/486*x^9*y^9*z^5 - 
      5/486*x^9*y^9*z^4 - 7/324*x^9*y^8*z^5 + 29/486*x^9*y^8*z^4 - 7/324*x^9*y^8*z^3 - 
      4/81*x^9*y^7*z^5 + 73/243*x^9*y^7*z^4 - 73/243*x^9*y^7*z^3 + 4/81*x^9*y^7*z^2 - 
      8/243*x^9*y^6*z^5 + 143/243*x^9*y^6*z^4 - 115/81*x^9*y^6*z^3 + 143/243*x^9*y^6*z^2 - 
      8/243*x^9*y^6*z + 5/486*x^9*y^5*z^5 - 7/54*x^9*y^5*z^4 + 25/243*x^9*y^5*z^3 - 25/243*x^9*y^5*z^2
      + 7/54*x^9*y^5*z - 5/486*x^9*y^5 - 145/486*x^9*y^4*z^4 + 722/243*x^9*y^4*z^3 - 
      1439/243*x^9*y^4*z^2 + 722/243*x^9*y^4*z - 145/486*x^9*y^4 - 109/972*x^9*y^3*z^3 + 
      565/972*x^9*y^3*z^2 - 565/972*x^9*y^3*z + 109/972*x^9*y^3 - 241/972*x^9*y^2*z^2 + 
      146/243*x^9*y^2*z - 241/972*x^9*y^2 - 23/243*x^9*y*z + 23/243*x^9*y - 8/243*x^9 - 
      1/486*x^8*y^10*z^5 - 5/486*x^8*y^9*z^5 + 5/486*x^8*y^9*z^4 + 31/1944*x^8*y^8*z^5 - 
      23/486*x^8*y^8*z^4 + 31/1944*x^8*y^8*z^3 - 31/1944*x^8*y^7*z^5 + 157/1944*x^8*y^7*z^4 - 
      157/1944*x^8*y^7*z^3 + 31/1944*x^8*y^7*z^2 - 53/972*x^8*y^6*z^5 + 617/972*x^8*y^6*z^4 - 
      655/486*x^8*y^6*z^3 + 617/972*x^8*y^6*z^2 - 53/972*x^8*y^6*z - 19/972*x^8*y^5*z^5 + 
      73/108*x^8*y^5*z^4 - 1553/486*x^8*y^5*z^3 + 1553/486*x^8*y^5*z^2 - 73/108*x^8*y^5*z + 
      19/972*x^8*y^5 + 115/972*x^8*y^4*z^4 - 1483/972*x^8*y^4*z^3 + 535/162*x^8*y^4*z^2 - 
      1483/972*x^8*y^4*z + 115/972*x^8*y^4 + 53/108*x^8*y^3*z^3 - 2347/972*x^8*y^3*z^2 + 
      2347/972*x^8*y^3*z - 53/108*x^8*y^3 + 23/216*x^8*y^2*z^2 - 65/243*x^8*y^2*z + 23/216*x^8*y^2 + 
      101/972*x^8*y*z - 101/972*x^8*y + 2/81*x^8 - 7/486*x^7*y^9*z^5 + 7/486*x^7*y^9*z^4 - 
      13/243*x^7*y^8*z^5 + 31/243*x^7*y^8*z^4 - 13/243*x^7*y^8*z^3 - 1/81*x^7*y^7*z^5 + 
      7/486*x^7*y^7*z^4 - 7/486*x^7*y^7*z^3 + 1/81*x^7*y^7*z^2 - 1/162*x^7*y^6*z^5 + 1/18*x^7*y^6*z^4 
      - 8/81*x^7*y^6*z^3 + 1/18*x^7*y^6*z^2 - 1/162*x^7*y^6*z - 7/486*x^7*y^5*z^5 + 
      199/486*x^7*y^5*z^4 - 428/243*x^7*y^5*z^3 + 428/243*x^7*y^5*z^2 - 199/486*x^7*y^5*z + 
      7/486*x^7*y^5 + 79/486*x^7*y^4*z^4 - 937/486*x^7*y^4*z^3 + 988/243*x^7*y^4*z^2 - 
      937/486*x^7*y^4*z + 79/486*x^7*y^4 - 83/162*x^7*y^3*z^3 + 451/162*x^7*y^3*z^2 - 
      451/162*x^7*y^3*z + 83/162*x^7*y^3 - 7/18*x^7*y^2*z^2 + 467/486*x^7*y^2*z - 7/18*x^7*y^2 - 
      11/243*x^7*y*z + 11/243*x^7*y - 5/243*x^7 - 73/1944*x^6*y^8*z^5 + 91/972*x^6*y^8*z^4 - 
      73/1944*x^6*y^8*z^3 - 61/648*x^6*y^7*z^5 + 295/648*x^6*y^7*z^4 - 295/648*x^6*y^7*z^3 + 
      61/648*x^6*y^7*z^2 - 1/36*x^6*y^6*z^5 + 151/648*x^6*y^6*z^4 - 407/972*x^6*y^6*z^3 + 
      151/648*x^6*y^6*z^2 - 1/36*x^6*y^6*z - 1/972*x^6*y^5*z^5 + 11/486*x^6*y^5*z^4 - 
      73/972*x^6*y^5*z^3 + 73/972*x^6*y^5*z^2 - 11/486*x^6*y^5*z + 1/972*x^6*y^5 + 29/324*x^6*y^4*z^4 
      - 103/108*x^6*y^4*z^3 + 35/18*x^6*y^4*z^2 - 103/108*x^6*y^4*z + 29/324*x^6*y^4 - 
      13/36*x^6*y^3*z^3 + 1889/972*x^6*y^3*z^2 - 1889/972*x^6*y^3*z + 13/36*x^6*y^3 + 
      655/972*x^6*y^2*z^2 - 1679/972*x^6*y^2*z + 655/972*x^6*y^2 + 149/972*x^6*y*z - 149/972*x^6*y + 
      7/972*x^6 - 43/972*x^5*y^7*z^5 + 73/324*x^5*y^7*z^4 - 73/324*x^5*y^7*z^3 + 43/972*x^5*y^7*z^2 - 
      31/486*x^5*y^6*z^5 + 641/972*x^5*y^6*z^4 - 643/486*x^5*y^6*z^3 + 641/972*x^5*y^6*z^2 - 
      31/486*x^5*y^6*z - 5/486*x^5*y^5*z^5 + 121/486*x^5*y^5*z^4 - 473/486*x^5*y^5*z^3 + 
      473/486*x^5*y^5*z^2 - 121/486*x^5*y^5*z + 5/486*x^5*y^5 + 1/243*x^5*y^4*z^4 - 37/972*x^5*y^4*z^3
      + 37/486*x^5*y^4*z^2 - 37/972*x^5*y^4*z + 1/243*x^5*y^4 - 46/243*x^5*y^3*z^3 + 
      233/243*x^5*y^3*z^2 - 233/243*x^5*y^3*z + 46/243*x^5*y^3 + 8/27*x^5*y^2*z^2 - 367/486*x^5*y^2*z 
      + 8/27*x^5*y^2 - 88/243*x^5*y*z + 88/243*x^5*y - 2/81*x^5 - 11/486*x^4*y^6*z^5 + 
      473/1944*x^4*y^6*z^4 - 161/324*x^4*y^6*z^3 + 473/1944*x^4*y^6*z^2 - 11/486*x^4*y^6*z - 
      7/486*x^4*y^5*z^5 + 383/972*x^4*y^5*z^4 - 1613/972*x^4*y^5*z^3 + 1613/972*x^4*y^5*z^2 - 
      383/972*x^4*y^5*z + 7/486*x^4*y^5 + 73/972*x^4*y^4*z^4 - 1499/1944*x^4*y^4*z^3 + 
      1525/972*x^4*y^4*z^2 - 1499/1944*x^4*y^4*z + 73/972*x^4*y^4 - 5/648*x^4*y^3*z^3 + 
      85/1944*x^4*y^3*z^2 - 85/1944*x^4*y^3*z + 5/648*x^4*y^3 + 175/972*x^4*y^2*z^2 - 
      145/324*x^4*y^2*z + 175/972*x^4*y^2 - 79/972*x^4*y*z + 79/972*x^4*y + 67/972*x^4 - 
      1/243*x^3*y^5*z^5 + 55/486*x^3*y^5*z^4 - 118/243*x^3*y^5*z^3 + 118/243*x^3*y^5*z^2 - 
      55/486*x^3*y^5*z + 1/243*x^3*y^5 + 13/162*x^3*y^4*z^4 - 143/162*x^3*y^4*z^3 + 149/81*x^3*y^4*z^2
      - 143/162*x^3*y^4*z + 13/162*x^3*y^4 - 49/243*x^3*y^3*z^3 + 509/486*x^3*y^3*z^2 - 
      509/486*x^3*y^3*z + 49/243*x^3*y^3 + 2/243*x^3*y^2*z^2 - 11/486*x^3*y^2*z + 2/243*x^3*y^2 - 
      20/243*x^3*y*z + 20/243*x^3*y + 1/54*x^2*y^4*z^4 - 43/216*x^2*y^4*z^3 + 5/12*x^2*y^4*z^2 - 
      43/216*x^2*y^4*z + 1/54*x^2*y^4 - 311/1944*x^2*y^3*z^3 + 559/648*x^2*y^3*z^2 - 559/648*x^2*y^3*z
      + 311/1944*x^2*y^3 + 479/1944*x^2*y^2*z^2 - 152/243*x^2*y^2*z + 479/1944*x^2*y^2 - 1/324*x^2*y*z
      + 1/324*x^2*y + 5/324*x^2 - 1/36*x*y^3*z^3 + 5/36*x*y^3*z^2 - 5/36*x*y^3*z + 1/36*x*y^3 + 
      43/324*x*y^2*z^2 - 55/162*x*y^2*z + 43/324*x*y^2 - 11/81*x*y*z + 11/81*x*y + 1/72*y^2*z^2 - 
      1/36*y^2*z + 1/72*y^2 - 1/27*y*z + 1/27*y + 1/36
  ;}

  {Lt=
  x^3*y^5*z^5*t^5 + 5/6*x^3*y^4*z^4*t^4 + 2/3*x^3*y^3*z^3*t^3 + 1/2*x^3*y^2*z^2*t^2 - 
      5/6*x^2*y^4*z^4*t^4 - 2/3*x^2*y^3*z^3*t^3 - 1/2*x^2*y^2*z^2*t^2 - 1/3*x^2*y*z*t + 
      2/3*x*y^3*z^3*t^3 + 1/2*x*y^2*z^2*t^2 + 1/3*x*y*z*t + 1/6*x - 1/2*y^2*z^2*t^2 - 1/3*y*z*t - 1/6
  ;}

  polrootsreal(res1y,[0,1])
  /* no roots */
.........................................................................

/* a = 1, g = 1:
  (1+x)*(1+x*y)*(1-x*y*z)*(1-2*x*y*z*t+(x*y*z*t)^2)*(1-x*y*z*t)*(1-y)*(1+y*z)*(1+2*y*z*t+(y*z*t)^2)*(1+y*z*t) */

  {res1y=
  x - 5/6
  ;}

  {res1z=
  x^2*y - 1/2*x*y + 5/6*x - 2/3
  ;}

  {res1t=
  x^3*y^2*z - 1/6*x^2*y^2*z + 5/6*x^2*y*z - 5/6*x^2*y - 1/3*x*y*z + 1/3*x*y - 2/3*x + 1/2
  ;}

  {Lt=
  x*y*z*t + 1/2*x - 1/2
  ;}

  polrootsreal(res1y,[0,1])
  x0 = 5/6
  
  polrootsreal(subst(res1z,x,x0),[0,1])
  /* no roots */
.........................................................................
