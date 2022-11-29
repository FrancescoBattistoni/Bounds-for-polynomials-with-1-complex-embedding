/* INNER */

  R(x,y,z,g) = (1+x)*(1-x*y)*(1+2*g*x*y+(x*y)^2)*(1+y)*(1-2*g*y*z+(y*z)^2);

  {resz=
  x^15 - 3109/2868*x^14 - 1505/2868*x^13 + 4493/11472*x^12 + 5675/5736*x^11 - 8887/11472*x^10 - 147/956*x^9 + 255/1912*x^8
  + 709/2868*x^7 - 1199/5736*x^6 + 73/2868*x^5 + 39/3824*x^4 + 107/5736*x^3 - 77/3824*x^2 + 7/956*x - 1/956
  ;}
  
  {res1y=
  x^5 + 3/2*x^4*z - 1/4*x^4 - 1/2*x^3*z - 1/2*x^3 - 1/2*x^2*z + 1/2*x - 1/4
  ;}
  
  {res1g=
  x^3*y^3 + 3/2*x^2*y^3*z + 3/4*x^2*y^3 - 3/4*x^2*y^2 + x*y^3*z - x*y^2*z - 1/2*x*y^2 + 1/2*x*y - 1/2*y^2*z + 1/4*y - 1/4
  ;}
  
  {Lg=
  x^2*y^2*z - x*y^2*z^2 + 4*x*y*z*g - x + z
  ;}
  
  v=polrootsreal(resz,[0,1]);  /* length 1 */
  x0=v[1];
  w=polrootsreal(subst(res1y,x,x0),[0,1]);  /* length 1 */
  z0=w[1];
  u=polrootsreal(subst(subst(res1g,x,x0),z,z0),[0,1]);  /* length 1 */
  y0=u[1]
  s=polrootsreal(subst(subst(subst(Lg,x,x0),z,z0),y,y0),[0,1]) /* length 1 */
  g0 = s[1]
  
  R(x0,y0,z0,g0)   /*  ===  2.0597880260425542589836374984271463532  */
....................................................

/* y = 1 */

  R(x,z,g) = (1+x)*(1-x)*(1+2*g*x+x^2)*2*(1-2*g*z+z^2);
  
  {resz=
  x^2 - 1/6
  ;}
  
  {res1g=
  x^3 + 3/2*x^2*z - 1/2*z
  ;}
  
  {Lx=
  x^3 + 3/2*x^2*g - 1/2*g
  ;}
  
  v=polrootsreal(resz,[0,1]);  /* length 1 */
  x0=v[1];
  w=polrootsreal(subst(res1g,x,x0),[0,1]);  /* length 1 */
  z0=w[1];
  u=polrootsreal(subst(subst(Lx,x,x0),z,z0),[0,1]);  /* length 1 */
  g0=u[1]
  
  R(x0,z0,g0)   /*  <=  2.1433470507544581618655692729766803842  */
....................................................
  
/* z = 0 */
  
  R(x,y)=(1+x)*(1-x*y)*(1+x*y)^2*(1+y);

  {resy=
  x^3 - 1/4*x^2 - 1/4
  ;}
  
  {Lx=
  x^2*y^2 + 3/4*x*y^2 - 1/4*x*y - 1/4*y - 1/4
  ;}
  
  v=polrootsreal(resy,[0,1]);  /* length 1 */
  x0=v[1];
  w=polrootsreal(subst(Lx,x,x0),[0,1]);
  y0=w[1];

  R(x0,y0)  /* <= 3.2854570820718519149100458947640051564  */
....................................................

/* z = 1 */

  {resy=
  x^3 - 262/125*x^2 + 33/25*x - 36/125
  ;}
  
  v=polrootsreal(resy,[0,1]);  /* no roots */
....................................................
  
/* z = 1, g = 0 */
  
  {resy=
  x^9 - 3/4*x^8 + 5/8*x^7 - 7/32*x^6 + 3/4*x^5 - 67/32*x^4 - 15/32*x^2 - 27/32
  ;}
  
  v=polrootsreal(resy,[0,1]);  /* no roots */
....................................................
 
/* z = 1, g = 1 */
  
  {resy=
  x^4 - 1/4*x^3 + 13/8*x^2 - 1/2*x + 3/8
  ;}
  
  v=polrootsreal(resy,[0,1]);  /* no roots */
....................................................
  
/* g = 0 */

  {resy=
  x^9 - 3/4*x^8 + 5/8*x^7 - 7/32*x^6 + 3/4*x^5 - 67/32*x^4 - 15/32*x^2 - 27/32
  ;}
  
  v=polrootsreal(resy,[0,1]);  /* no roots */
....................................................
  
/* g = 1 */
  
  R(x,y) = (1+x)*(1-x*y)*(1+x*y)^2*(1+y);
  
  {resy=
  x^3 - 1/4*x^2 - 1/4
  ;}

  {Lx=
  x^2*y^2 + 3/4*x*y^2 - 1/4*x*y - 1/4*y - 1/4
  ;}
  
  v=polrootsreal(resy,[0,1]);  /* length 1 */
  x0=v[1];
  w=polrootsreal(subst(Lx,x,x0),[0,1]);
  y0=w[1];
  
  R(x0,y0)  /* <= 3.2854570820718519149100458947640051564  */
....................................................
  
