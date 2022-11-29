

  {res1y=
  x - 1/2
  ;}
  
  {res1z=
  x*y + 1/2*y - 1/2
  ;}
  
  {Lz=
  x*y*z - 1/2
  ;}
  
  v=polrootsreal(res1y,[0,1]) /* length 1 */
  x0=v[1]

  u= polrootsreal(subst(res1z,x,1/2),[0,1]) /* length 1 */
  y0=u[1] 
  w= polrootsreal(subst(subst(Lz,x,x0),y,y0),[0,1])
  /* no roots in the open set */
....................................................

/* z= 1: (1+x)*(1-x*y)*(1-x*y)*(1+x*y)^3*(1+y) */
  
  R(x,y)=(1+x)*(1-x*y)*(1-x*y)*(1+x*y)^3*(1+y)
  
  {res1y=
   x^3 - 1/6*x^2 - 1/6
  ;}  

  {Ly=
  x^2*y^2 + 5/6*x^2*y - 1/6*x*y - 1/6*x - 1/6
  ;}

  v=polrootsreal(res1y,[0,1])  /* length 1 */
  x0=v[1]
  u=polrootsreal(subst(Ly,x,x0),[0,1])
  y0=u[1]
  R(x0,y0)
  \\ =  2.6399736753739488164926638942353861445
....................................................
