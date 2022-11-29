/* x = 1 */

  y0=sqrt(4/17);

  Ly=y^3 + 3/2*y^2*g - 1/2*g;

  R(y,g)=2*(1-y)*(1+2*y*g+y^2)*(1+y)*(1-g^2)^(1/8);

  v=polrootsreal(subst(Ly,y,y0),[0,1]);
  g0=v[1];

  R(y0,g0)
  \\ = 2.7097512625653105400215894723018612937
....................................................

