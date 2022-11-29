/* INTERIOR */

  L=(1+t^2)*(1+x*y*t)*(1-x*y*t*a)*(1+y^2)*(1-y*t)*(1+y*t*a);

  {res1y=
  x^4 - x^3 + 7/4*x^2 - 7/4*x + 1/2;
  }

  {res1t=
  x^2*y^2 - x^2 + 3*x*y^2 + 2*x - 2*y^2 - 1;
  }

  {res1a=
  x*y^2*t^2 - 2/3*x*y*t + 1/3*x*t^2 + 2/3*y^2*t^2 + 1/3*t^2 - 1/3;
  }

  {Lx=
  x*y*t*a + 1/2*a - 1/2;
  }

  v = polrootsreal(res1y,[0,1])  /* length 2 */
  x0=v[1];
  u = polrootsreal(subst(res1t,x,x0),[0,1])  /* no roots */

  x0=v[2];
  u = polrootsreal(subst(res1t,x,x0),[0,1])  /* length 1 */
  y0=u[1];
  w = polrootsreal(subst(subst(res1a,x,x0),y,y0),[0,1])  /* length 1 */
  t0=w[1];
  ua = polrootsreal(subst(subst(subst(Lx,x,x0),y,y0),t,t0),[0,1])  /* length 1 */ 
  a0=ua[1];
  subst(subst(subst(subst(L,x,x0),y,y0),t,t0),a,a0)
  \\ 1.5652995727446142266427366734658525993       <--------  value 
.........................................................................
  
