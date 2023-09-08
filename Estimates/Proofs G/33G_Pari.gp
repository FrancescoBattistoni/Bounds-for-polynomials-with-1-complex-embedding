
/* x = 0, y = 1 */
  
  R(z,g) = 2*(1+g)*(1-z)*(1+z)*(1+2*g*z+z^2)*2*sqrt(1-g^2)

  {res1g=
  z^3 + 1/4*z^2 - 1/6*z - 1/12
  ;}
  
  {Lz=
  z^3 + 3/2*z^2*g - 1/2*g
  ;}
  
  v=polrootsreal(res1g,[0,1]);  /* length one */
  z0=v[1];
  u=polrootsreal(subst(Lz,z,z0),[0,1]);  /* length one */
  g0=u[1];
  R(z0,g0)  \\    7.1586113968778701267626122275096609752
....................................................
  
/* g = 0 */
  
  L= (1+x*y*z)*(1-x*y*z*t)*(1+y^2)*(1-y*z)*(1+y*z*t)*(1+z^2)*2
  
  {resy=
  x^7 + 2*x^6 + 9/4*x^5 + 9/4*x^4 - 3/4*x^3 - 3/4*x^2 - 2*x + 1
  ;}
  
  {res1z=
  x^4*y^2 - 1/4*x^2 - 3/4*x + 1/2
  ;}

  {res1t=
  x^2*y*z + 1/2*x - 1/2
  ;}

  {Lx=
  x*y*z*t + 1/2*t - 1/2
  ;}

  vx=polrootsreal(resy,[0,1]);  /* length two */
  x0=vx[1];
  vy=polrootsreal(subst(res1z,x,x0),[0,1]);  /* no roots */
  
  x0=vx[2];
  vy=polrootsreal(subst(res1z,x,x0),[0,1]);  /* length one */
  y0=vy[1];
  vz=polrootsreal(subst(subst(res1t,y,y0),x,x0),[0,1]);  /* length one */
  z0=vz[1];
  vt=polrootsreal(subst(subst(subst(Lx,z,z0),y,y0),x,x0),[0,1]);  /* length one */
  t0=vt[1];
  my(x=x0,y=y0,z=z0,t=t0);eval(L)   \\   3.1305991454892284532854733469317051985
....................................................
 
/* y = 1 */

   R(x,z,t,g)= (1+x*z)*(1-x*z*t)*2*(1+g)*(1-z)*(1+z*t)*(1+2*g*z+z^2)*2*sqrt(1-g^2)
  
  {resz=
  x^5 - 11/6*x^4 + 7/3*x^3 - 3*x^2 + 17/6*x - 1
  ;}
  
  {res1t=
  x^2*z + 1/2*x - 1/2
  ;}
  
  {res1g=
  x*z*t + 1/2*t - 1/2
  ;}
  
  {Lg=
  z^2*g - 1/2*z^2 + 3*z*g^2 - z*g - z + g - 1/2
  ;}

  vx=polrootsreal(resz,[0,1]);  /* length one */
  x0=vx[1];
  vz=polrootsreal(subst(res1t,x,x0),[0,1]);  /* length one */
  z0=vz[1];
  vt=polrootsreal(subst(subst(res1g,z,z0),x,x0),[0,1]);  /* length one */
  t0=vt[1];
  vg=polrootsreal(subst(subst(subst(Lg,z,z0),y,y0),x,x0),[0,1]);  /* length one */
  g0=vg[1];
  R(x0,z0,t0,g0)   \\      R(x,z,t,g)= (1+x*z)*(1-x*z*t)*2*(1+g)*(1-z)*(1+z*t)*(1+2*g*z+z^2)*2*sqrt(1-g^2)
  
  {resz=
  x^5 - 11/6*x^4 + 7/3*x^3 - 3*x^2 + 17/6*x - 1
  ;}
  
  {res1t=
  x^2*z + 1/2*x - 1/2
  ;}
  
  {res1g=
  x*z*t + 1/2*t - 1/2
  ;}
  
  {Lg=
  z^2*g - 1/2*z^2 + 3*z*g^2 - z*g - z + g - 1/2
  ;}

  vx=polrootsreal(resz,[0,1]);  /* length one */
  x0=vx[1];
  vz=polrootsreal(subst(res1t,x,x0),[0,1]);  /* length one */
  z0=vz[1];
  vt=polrootsreal(subst(subst(res1g,z,z0),x,x0),[0,1]);  /* length one */
  t0=vt[1];
  vg=polrootsreal(subst(subst(subst(Lg,z,z0),y,y0),x,x0),[0,1]);  /* length one */
  g0=vg[1];
  R(x0,z0,t0,g0)   \\   6.6800672476067084452776060856129284315
....................................................
  
/* g = 1 */

  R(x,y,z,t,g) = (1+x*y*z)*(1-x*y*z*t)*(1+2*g*y+y^2)*(1-y*z)*(1+y*z*t)*(1+2*g*z+z^2)*2*sqrt(1-g^2)
  
  {resy=
  x^10 - 6*x^9 + 671/36*x^8 - 1567/36*x^7 + 4223/72*x^6 - 469/8*x^5 + 793/72*x^4 + 5015/72*x^3 -
    3029/36*x^2 + 113/3*x - 6
  ;}

  {res1z=
  x^7*y^4 - 16/13*x^7*y^2 + 3/13*x^7 + 54/13*x^6*y^4 + 77/26*x^6*y^2 - 35/26*x^6 - 15/13*x^5*y^4 -
    171/52*x^5*y^2 + 171/52*x^5 - 120/13*x^4*y^4 + 115/26*x^4*y^2 - 227/52*x^4 + 8*x^3*y^4 -
    183/26*x^3*y^2 + 89/26*x^3 - 24/13*x^2*y^4 + 177/26*x^2*y^2 - 21/13*x^2 - 167/52*x*y^2 + 23/52*x
    + 15/26*y^2 - 3/52
  ;}

  {res1t=
  x^2*y*z + 1/2*x - 1/2
  ;}

  {res4g=
  x*y*z*t + 1/2*x - 1/2
  ;}

  {Lg=
  y^2*z^2*g + 4*y^2*z*g^2 - 2*y^2*z + y^2*g + 4*y*z^2*g^2 - 2*y*z^2 + 12*y*z*g^3 - 8*y*z*g + 4*y*g^2 -
    2*y + z^2*g + 4*z*g^2 - 2*z + g
  ;}

  vx=polrootsreal(resy,[0,1]);  /* length two */
  x0=vx[1]; \\ this is 1/2
  vy=polrootsreal(subst(res1z,x,x0),[0,1]);  /* length one but the root is 1 so it is not in the open set */

  x0=vx[2];
  vy=polrootsreal(subst(res1z,x,x0),[0,1]);  /* length one */
  y0=vy[1];
  vz=polrootsreal(subst(subst(res1t,y,y0),x,x0),[0,1]);  /* length one */
  z0=vz[1];
  vt=polrootsreal(subst(subst(subst(res4g,z,z0),y,y0),x,x0),[0,1]);  /* length one */
  t0=vt[1];
  vg=polrootsreal(subst(subst(subst(subst(Lg,t,t0),z,z0),y,y0),x,x0),[0,1]);  /* length one */
  g0=vg[1];
  R(x0,y0,z0,t0,g0) \\  6.1686245360610869196334305419235384352
....................................................
