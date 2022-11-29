/* x = 0, t = 1 */

  {resz=
  y^11 - 1/3*y^10 + 9/17*y^9 - 3/17*y^8 - 425542/585225*y^7 + 425542/1755675*y^6 -
      695446/2926125*y^5 + 695446/8778375*y^4 + 10949/65025*y^3 - 10949/195075*y^2
      - 13/375*y + 13/1125
  ;}
  
  {res1g=
  y^4*z^3 - y^3*z^2 - 20/27*y^2*z^5 - 2/9*y^2*z^3 + 5/27*y^2*z + 16/27*y*z^4 + 
      10/27*y*z^2 + 1/27*y - 2/27*z^3 - 4/27*z
  ;}
  
  {Ly=
  y^2*z + 4/3*y*z*g - 2/3*y + 1/3*z - 2/3*g
  ;}
  
  R(y,z,g)=(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*(1-2*z*g+z^2)*2*(1-g^2)^(1/4); 
  
  v=polrootsreal(resz,[0,1]);  /* Length 2 */
  y0=v[1];
  w=polrootsreal(subst(res1g,y,y0),[0,1]);  /* Length 1, since the second root is z = 1 */
  z0=w[1];
  u=polrootsreal(subst(subst(Ly,y,y0),z,z0),[0,1]);  /* No roots */
  
  
  y0=v[2];
  w=polrootsreal(subst(res1g,y,y0),[0,1]);  /* Length 3 */
  z0=w[1];
  u=polrootsreal(subst(subst(Ly,y,y0),z,z0),[0,1]);  /* No roots */

  z0=w[2];
  u=polrootsreal(subst(subst(Ly,y,y0),z,z0),[0,1]);  /* Length 1 */
  g0=u[1];
  R(y0,z0,g0) /* < 3.7 */
  
  z0=w[3];
  u=polrootsreal(subst(subst(Ly,y,y0),z,z0),[0,1]);  /* No roots */
.........................................................................
  
/* x = 1, t = 0 */

  {resz=
  y^11 - 1/2*y^10 + 393/425*y^9 - 393/850*y^8 + 63838/195075*y^7 - 
      31919/195075*y^6 - 434/39015*y^5 + 217/39015*y^4 - 313/7803*y^3 + 
      313/15606*y^2 - 65/7803*y + 65/15606
  ;}
  
  {res1g=
  y^5*z^2 - 4/5*y^4*z - 27/20*y^3*z^4 + 3/10*y^3*z^2 + 1/10*y^3 + 27/20*y^2*z^3 - 
      1/2*y^2*z - 1/4*y*z^2 + 1/5*y - 1/20*z
  ;}
  
  {Lz=
  y*z^2 + 4/3*y*z*g + 1/3*y - 2/3*z - 2/3*g
  ;}
  
  R(y,z,g)=(1-2*y*g+y^2)*(1+2*y*g+y^2)*(1-y*z)*(1+2*z*g+z^2)*2*(1-g^2)^(1/4);
  
  v=polrootsreal(resz,[0,1]);  /* Length 2 */
  y0=v[1];
  w=polrootsreal(subst(res1g,y,y0),[0,1]);  /* Length 2 */
  z0=w[1];
  u=polrootsreal(subst(subst(Lz,y,y0),z,z0),[0,1]);  /* No roots */

  z0=w[2];
  u=polrootsreal(subst(subst(Lz,y,y0),z,z0),[0,1]);  /* No roots */
  
  y0=v[2];
  w=polrootsreal(subst(res1g,y,y0),[0,1]);  /* Length 1 */
  z0=w[1];
  u=polrootsreal(subst(subst(Lz,y,y0),z,z0),[0,1]);  /* Length 1 */
  g0=u[1];
  
  R(y0,z0,g0)  /* < 3.7 */
.........................................................................
  
  
  
