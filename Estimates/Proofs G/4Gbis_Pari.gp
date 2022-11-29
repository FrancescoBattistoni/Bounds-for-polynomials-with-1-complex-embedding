/* Study of inequality 4Gbis

      - - +                 <= 4.36 for y in [0,1/10]
          -
          -

*/

  L=(1+x)*(1+x*y)*(1-x*y*z)*(1+y*z)*(1+z);


/* left-boundaries */

/* x = 0: (1+y*z)*(1+z). Maximized at y = z = 1 giving the value 4 */

/* y = 0: (1+x)*(1+z).   Maximized at x = z = 1 giving the value 4 */

/* z = 0: (1+x)*(1+x*y). Maximized at x = y = 1 giving the value 4 */


/* right-boundaries */

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
  v=polrootsreal(subst(Lz,z,z0),[0,1/10]) /* no roots */

  /* y = 0: already discussed in greater detail (it is <= 4) */
  /* z = 0: already discussed in greater detail (it is <= 4) */
  /* y = 1/10: 2*(1+y)*(1-(y*z)^2)*(1+z) <= 4.36 */ /* PARI my(y=1/10);ploth(z=0,1   ,2*(1+y)*(1-(y*z)^2)*(1+z)) */
  /* z = 1   : 2*(1+y)*(1-(y*z)^2)*(1+z) <= 4.36 */ /* PARI my(z=1)   ;ploth(y=0,1/10,2*(1+y)*(1-(y*z)^2)*(1+z)) */


/* y = 1/10: (1+x)*(1+x*(1/10))*(1-x*(1/10)*z)*(1+(1/10)*z)*(1+z) */
  L=(1+x)*(1+x*(1/10))*(1-x*(1/10)*z)*(1+(1/10)*z)*(1+z);
  Tx=factor(deriv(L,x))
  Tz=factor(deriv(L,z))
  Lx=Tx[3,1];
  Lz=Tz[3,1];

  T=factor(polresultant(Lx,Lz,x))
  res1x=T[3,1]

  u=polrootsreal(res1x,[0,1]) /* no roots */

  /* x = 0: already discussed in greater detail (it is <= 4) */
  /* z = 0: already discussed in greater detail (it is <= 4) */
  /* x = 1: (1+x)*(1+x*(1/10))*(1-x*(1/10)*z)*(1+(1/10)*z)*(1+z) <= 4.36 */ /* PARI my(x=1);ploth(z=0,1,(1+x)*(1+x*(1/10))*(1-x*(1/10)*z)*(1+(1/10)*z)*(1+z)) */
  /* z = 1: (1+x)*(1+x*(1/10))*(1-x*(1/10)*z)*(1+(1/10)*z)*(1+z) <= 4.36 */ /* PARI my(z=1);ploth(x=0,1,(1+x)*(1+x*(1/10))*(1-x*(1/10)*z)*(1+(1/10)*z)*(1+z)) */


/* z = 1: (1+x)*(1+x*y)*(1-x*y)*(1+y)*2 */
  L=(1+x)*(1+x*y)*(1-x*y)*(1+y)*2;
  Tx=factor(deriv(L,x))
  Ty=factor(deriv(L,y))
  Lx=Tx[2,1];
  Ly=Ty[2,1];

  T=factor(polresultant(Lx,Ly,x))
  res1x=T[3,1];

  u=polrootsreal(res1x,[0,1/10]) /* no roots */

  /* x = 0: already discussed in greater detail (it is <= 4) */
  /* y = 0: already discussed in greater detail (it is <= 4) */
  /* x = 1: already discussed in greater detail (it is <= 4.36) */
  /* y = 1: already discussed in greater detail (it is <= 4.36) */

/* The function is <=  4.36 */
