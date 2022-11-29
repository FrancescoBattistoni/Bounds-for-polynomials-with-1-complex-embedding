/* a = 1, Interior */

  {rest=
  z^13 - 1/3*z^12 + 14/405*z^11 - 14/1215*z^10 - 55907/54675*z^9 + 
      55907/164025*z^8 + 41996/273375*z^7 - 41996/820125*z^6 + 378133/1366875*z^5 
      - 378133/4100625*z^4 - 1802/18225*z^3 + 1802/54675*z^2 + 1/50625*z - 
      1/151875
  ;}
  
  {res1y=
  z^4*t^4 - 2/3*z^3*t^5 - 2/3*z^3*t^3 - 4/3*z^2*t^6 + 5/9*z^2*t^4 + 1/9*z^2*t^2 + 
      10/9*z*t^5 + 4/9*z*t^3 - 2/9*z*t - 1/9*t^4 - 1/3*t^2 + 1/9
  ;}
  
  {res1g=
  y*z^2*t^2 - 1/2*y*z*t + 1/2*z*t + 1/3*t^2 - 1/3
  ;}
  
  {Ly=
  y^2*z^2*t - 4/3*y*z*t*g + 2/3*y*z + 1/3*t - 2/3*g
  ;}
  
  R(y,z,t,g)=(1+y*z*t)*(1-2*y*z*g+(y*z)^2)*(1+2*z*g+z^2)*(1-z*t)*(1+2*t*g+t^2)*(1-2*t*g+t^2)*2*sqrt(1-g^2);
  
  v=polrootsreal(rest,[0,1]);  /* Length 3 */
  z0=v[1];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* No roots */
  
  z0=v[2];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* No roots */
  
  z0=v[3];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 2 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* Length 1 */
  y0=u[1];
  ug=polrootsreal(subst(subst(subst(Ly,t,t0),z,z0),y,y0),[0,1]);  /* Length 1 */
  g0=ug[1];
  R(y0,z0,t0,g0)  /* < 3.61 */
  
  t0=w[2];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* Length 1 */
  y0=u[1];
  ug=polrootsreal(subst(subst(subst(Ly,t,t0),z,z0),y,y0),[0,1]);  /* Length 1 */
  g0=ug[1];
  R(y0,z0,t0,g0)  /* < 3.33 */
....................................................  
  
/* a = 0, Interior */
  
  {rest=
  z^4 - 22/35*z^2 + 1/25
  ;}
  
  {res1y=
  z^2*t^2 - 8/3*z*t^3 + 2/3*z*t + 5/3*t^2 - 2/3
  ;}
  
  {res1g=
  y*z^2*t^2 - 1/2*y*z*t + 1/2*z*t + 1/3*t^2 - 1/3
  ;}

  {Ly=
  y^2*z^2*t - 4/3*y*z*t*g + 2/3*y*z + 1/3*t - 2/3*g
  ;}
  
  R(y,z,t,g)=(1+y*z*t)*(1-2*y*z*g+(y*z)^2)*(1+2*z*g+z^2)*(1-z*t)*(1+2*t*g+t^2)*2*sqrt(1-g^2);
  v=polrootsreal(rest,[0,1]);  /* Length 2 */
  z0=v[1];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* No roots */
  
  z0=v[2];
  w=polrootsreal(subst(res1y,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(res1g,t,t0),z,z0),[0,1]);  /* Length 1 */
  y0=u[1];
  ug=polrootsreal(subst(subst(subst(Ly,t,t0),z,z0),y,y0),[0,1]);  /* Length 1 */
  g0=ug[1];
  R(y0,z0,t0,g0)  /* < 4.1 */
....................................................  
  
/* a = 0, y = 1  */
  
  {rest=
  z^10 + 13/33*z^8 - 86/825*z^6 - 2/15*z^4 - 1/275*z^2 + 3/275
  ;}
  
  {res1g=
  z^6*t^4 - 32/27*z^4*t^6 + 4/27*z^4*t^4 - 17/27*z^4*t^2 + 25/27*z^2*t^4 - 
      2/9*z^2*t^2 + 2/27*z^2 - 5/27*t^2 + 2/27
  ;}
  
  {Lt=
  z^2*t^3 + 3/2*z^2*t^2*g + 1/2*z^2*t - 1/2*t - 1/2*g
  ;}

  R(z,t,g)=(1+z*t)*(1-2*z*g+z^2)*(1+2*z*g+z^2)*(1-z*t)*(1+2*t*g+t^2)*2*sqrt(1-g^2);
  
  v=polrootsreal(rest,[0,1]);  /* Length 2 */
  z0=v[1];
  w=polrootsreal(subst(res1g,z,z0),[0,1]);  /* Length 1 */
  t0=w[1];
  u=polrootsreal(subst(subst(Lt,t,t0),z,z0),[0,1]);  /* No roots */
  
  z0=v[2];
  w=polrootsreal(subst(res1g,z,z0),[0,1]);  /* Length 3 */
  t0=w[1];
  u=polrootsreal(subst(subst(Lt,t,t0),z,z0),[0,1]);  /* No roots */
  t0=w[2];
  u=polrootsreal(subst(subst(Lt,t,t0),z,z0),[0,1]);  /* No roots */
  t0=w[3];
  u=polrootsreal(subst(subst(Lt,t,t0),z,z0),[0,1]);  /* Length 1 */
  g0=u[1];
  R(z0,t0,g0)  /* < 5.5 */
....................................................
