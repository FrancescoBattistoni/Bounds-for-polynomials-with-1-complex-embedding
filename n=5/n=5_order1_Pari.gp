
\\ CASE 3 

{resy=
x^12 + 25/14*x^11 - 5433/2744*x^10 - 353053/27440*x^9 - 15563347/614656*x^8 - 
    24815439/878080*x^7 - 112139957/6146560*x^6 - 6710161/1536640*x^5 + 
    2526387/768320*x^4 + 21728727/6146560*x^3 + 8724483/6146560*x^2 + 
    173583/614656*x + 4131/153664
;}

{res1zc=
x^3*y^4 + 1/2*x^3*y^3 - 13/4*x^2*y^4 + 3/8*x^2*y^3 - 3/4*x^2*y^2 - 7/8*x^2*y + 
    3/2*x*y^4 + x*y^3 + 3/2*x*y^2 + 1/2*x - 3/8*y^3 - 3/4*y^2 - 5/8*y + 1/4
;}

{res1g=
x^3*y^4*z^2 - 13/4*x^2*y^4*z^2 + x^2*y^3*z^2 - 3/4*x^2*y + 3/2*x*y^4*z^2 + 
    5/4*x*y^3*z^2 - 5/4*x*y^2*z^2 + 5/4*x*y^2 - 1/4*x*y + 1/2*x - 3/4*y^3*z^2 + 
    1/2*y^2*z^2 - 1/2*y^2 - 1/2*y + 1/4
;}

{Lx=
x^3*y^3*z^2 - 3/4*x^2*y^3*z^2 - 3/4*x^2*y^2*z^2 - 3/2*x^2*y^2*z*g + 
    1/2*x*y^2*z^2 + x*y^2*z*g + x*y*z*g + 1/2*x*y - 1/2*y*z*g - 1/4*y - 1/4
;}

R(x,y,z,g)=(1-x)*(1-x*y)*(1-2*x*y*z*g+(x*y*z)^2)*(1-y)*(1-2*y*z*g+(y*z)^2)*(1-2*z*g+z^2)*2*sqrt(1-g^2);

/* We use the following program */

search(vx)={
  for(i=1,#vx,
     vy = polrootsreal(subst(res1zc,x,vx[i]),[-1,1]);
     if (#vy > 0,
        vz = vector(#vt,r,100);
        vg = vector(#vt,r,100);
	vR = vector(#vt);

        print("x= ",vx[i]);
        print("y z g R(x,y,z,g)");

        for(j=1,#vy,

              vz[j] = polrootsreal(subst(subst(res1g,x,vx[i]),y,vy[j]),[-1,1]);

              if(#vz[j]>0,
                 a = vector(#vz[j]);
                 b = vector(#vz[j]);

                 for(k=1,#vz[j], a[k]=polrootsreal(subst(subst(subst(Lx,x,vx[i]),y,vy[j]),z,vz[j][k]))[1];

                     if(abs(a[k])< 1, b[k]=R(vx[i],vy[j],vz[j][k],a[k]); 
                       )
                     ); 
                 vg[j]=a;
                 vR[j]=b;
                 );
              
          print("\n y=",vy[j],"\n z=", vz[j],"\n g=",vg[j],"\n R=",vR[j],"\n");
            );
         );
      );
}

/* Define v = polrootsreal(resy,[-1,1]) and apply search(v). You will obtain the result (i.e. the points we found give values under 16*M ) unless for the second root, which is -1/2: there you have to perform hand subtitutions, but you will find z=0 */
\\ END CASE 3
