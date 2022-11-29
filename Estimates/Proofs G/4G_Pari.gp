/* Study of inequality 4G

      - - +                 <= 128/27
          -
          -

*/

  L=(1+x)*(1+x*y)*(1-x*y*z)*(1+y*z)*(1+z);


/* 0-boundaries */

/* x = 0: (1+y*z)*(1+z). Maximized at y=z=1 giving the value 4 */

/* y = 0: (1+x)*(1+z). Maximized at x=z=1 giving the value 4 */

/* z = 0: (1+x)*(1+x*y). Maximized at x=y=1 giving the value 4 */


/* 1-boundaries */

/* x = 1: 2*(1+y)*(1-(y*z)^2)*(1+z)  */
  L=2*(1+y)*(1-(y*z)^2)*(1+z);
  Ty=factor(deriv(L,y))
  Tz=factor(deriv(L,z))
  Ly=Ty[2,1];
  Lz=Tz[2,1];

  T=factor(polresultant(Ly,Lz,y))
  resy=T[3,1];

  u=polrootsreal(resy,[0,1]) /* length 1 */
  z0=u[1]
  v=polrootsreal(subst(Lz,z,z0),[0,1]) /* length 1 */
  y0=v[1]
  subst(subst(L,y,y0),z,z0)
  \\ 4.4771685631181802012741118480453117551  <----

  /* y = 0: already discussed in greater detail (it is <= 4) */
  /* z = 0: already discussed in greater detail (it is <= 4) */
  /* y = 1: 2*(1+y)*(1-(y*z)^2)*(1+z) = 2*2*(1-z^2)*(1+z) <= 128/27 */
  /* z = 1: 2*(1+y)*(1-(y*z)^2)*(1+z) = 2*(1+y)*(1-y^2)*2 <= 128/27 */


/* y = 1: (1+x)^2*(1-x*z)*(1+z)^2  */
  L=(1+x)^2*(1-x*z)*(1+z)^2;
  Tx=factor(deriv(L,x))
  Tz=factor(deriv(L,z))
  Lx=Tx[3,1];
  Lz=Tz[3,1];

  T=factor(polresultant(Lx,Lz,x));
  res1x=T[2,1]

  u=polrootsreal(res1x,[0,1]) /* length 1 */
  z0=u[1]
  v=polrootsreal(subst(Lz,z,z0),[0,1]) /* length 1 */
  x0=v[1]
  subst(subst(L,x,x0),z,z0)
  \\ 4.2866941015089163237311385459533607682  <----

  /* x = 0: already discussed in greater detail (it is <= 4) */
  /* z = 0: already discussed in greater detail (it is <= 4) */
  /* x = 1: (1+x)^2*(1-x*z)*(1+z)^2 <= 128/27 */ /* PARI my(x=1);ploth(z=0,1,(1+x)^2*(1-x*z)*(1+z)^2) */
  /* z = 1: (1+x)^2*(1-x*z)*(1+z)^2 <= 128/27 */ /* PARI my(z=1);ploth(x=0,1,(1+x)^2*(1-x*z)*(1+z)^2) */


/* z = 1: (1+x)*(1+x*y)*(1-x*y)*(1+y)*2 */
  L=(1+x)*(1+x*y)*(1-x*y)*(1+y)*2;
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Lx=Tx[2,1];
  Ly=Ty[2,1];

  T=factor(polresultant(Lx,Ly,x))
  res1x=T[3,1];

  u=polrootsreal(res1x,[0,1]) /* length 1 */
  y0=u[1]
  v=polrootsreal(subst(Lx,y,y0),[0,1]) /* length 1 */
  x0=v[1]
  subst(subst(L,x,x0),y,y0)
  \\ 4.4771685631181802012741118480453117551  <----

  /* x = 0: already discussed in greater detail (it is <= 4) */
  /* y = 0: already discussed in greater detail (it is <= 4) */
  /* x = 1: already discussed in greater detail (it is <= 128/27) */
  /* y = 1: already discussed in greater detail (it is <= 128/27) */

/* The function is then <=  128/27 */
