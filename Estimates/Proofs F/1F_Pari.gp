/* z = 1 */


{rest=
y^5 - 487/648*y^4 + 5/36*y^3 - 38/81*y^2 + 29/108*y - 25/648
;}


{res1g=
y^2*t^2 - 2/3*y^2*t + 8/9*y*t^3 - 2/3*y*t^2 + 7/9*y*t - 1/3*y + 5/9*t^2 - 1/3*t 
    + 1/9
;}


{Ly=
y^2*t - 4/3*y*t*g + 2/3*y + 1/3*t - 2/3*g
;}


R(y,t,g)=4*(1-2*y*g+y^2)*(1+y*t)*(1+g)*(1-t)*(1+2*t*g+t^2)*sqrt(1-g^2);

v=polrootsreal(rest,[0,1]);  /* length 1 */

y0=v[1];

w=polrootsreal(subst(res1g,y,y0),[0,1]);

t0=w[1];

u=polrootsreal(subst(subst(Ly,t,t0),y,y0));

g0=u[1];

R(y0,t0,g0)   /* < 0.942 */


/* y = 1 */

{rest=
z^8 + 8/11*z^6 + 38/275*z^4 - 24/275*z^2 - 9/275
;}


{res1g=
z^6*t^4 - 32/27*z^4*t^6 + 4/27*z^4*t^4 - 17/27*z^4*t^2 + 25/27*z^2*t^4 - 
    2/9*z^2*t^2 + 2/27*z^2 - 5/27*t^2 + 2/27
;}


{Lt=
z^2*t^3 + 3/2*z^2*t^2*g + 1/2*z^2*t - 1/2*t - 1/2*g
;}

R(z,t,g)=2*(1-2*z*g+z^2)*(1+z*t)*(1+2*z*g+z^2)*(1-z*t)*(1+2*t*g+t^2)*sqrt(1-g^2);

v=polrootsreal(rest,[0,1]);  /* length 1 */

z0=v[1];

w=polrootsreal(subst(res1g,z,z0),[0,1]);  /* length 3 */

t0=w[1];

u=polrootsreal(subst(subst(Lt,t,t0),z,z0));

g0=u[1];  /* This is negative */


t0=w[2];

u=polrootsreal(subst(subst(Lt,t,t0),z,z0));

g0=u[1];  /* This is negative */


t0=w[3];

u=polrootsreal(subst(subst(Lt,t,t0),z,z0));

g0=u[1];

R(z0,t0,g0)   /* < 5.48 */


/* Interior */

R(y,z,t,g)=(1-2*y*z*g+(y*z)^2)*(1+y*z*t)*(1+2*z*g+z^2)*(1-z*t)*(1+2*t*g+t^2)*2*sqrt(1-g^2);

{resy=
z^4 - 22/35*z^2 + 1/25
;}


{res1t=
y^3*z^4 - 8/3*y^2*z^4 + 5/3*y^2*z^2 + y*z^4 - 2*y*z^2 + 5/3*z^2 - 2/3
;}


{res1g=
y*z^2*t^2 - 1/2*y*z*t + 1/2*z*t + 1/3*t^2 - 1/3
;}


{Ly=
y^2*z^2*t - 4/3*y*z*t*g + 2/3*y*z + 1/3*t - 2/3*g
;}

v=polrootsreal(resy,[0,1]);  /* length 2 */


z0=v[1];

w=polrootsreal(subst(res1t,z,z0),[0,1]);  /* No roots */

z0=v[2];

w=polrootsreal(subst(res1t,z,z0),[0,1]);  /* Length 1 */

y0=w[1];

u=polrootsreal(subst(subst(res1g,y,y0),z,z0),[0,1]);

t0=u[1];

ug=polrootsreal(subst(subst(subst(Ly,y,y0),z,z0),t,t0),[0,1]);

g=g0;

R(y0,z0,t0,g0)  /* < 3.82 */



