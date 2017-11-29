function [xmin,fxmin]=brent(ax,bx,cx,tol,xfunc,yfunc)
    %%%% Written by Gaoyang Li, 23/Nov./2017; 
    %%%% Adapted ("translated") from Numerical Recipes for C;
    
  iter=100;
  cgold=0.3819660;
  zeps=1e-18;

  a=min(ax,cx);
  b=max(ax,cx);
  
  v=bx;
  w=v;
  x=v;
  e=0.0;
  fx=ls(x,xfunc,yfunc);
  fv=fx;
  fw=fx;
  
  for i=1:iter;
      xm=0.5*(a+b);
      tol1=tol*abs(x)+zeps;
      tol2=2.0*tol1;
      
      if (abs(x-xm) <= (tol2-0.5*(b-a)));
         break;    
      end;
      if (abs(e) >= tol1);
          r=(x-w)*(fx-fv);
          q=(x-v)*(fx-fw);
          p=(x-v)*q-(x-w)*r;
          q=2.0*(q-r);
          if (q>=0);
              p=-p;
          end;
          q=abs(q);
          etemp=e;
          e=d;
          if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
              if (x>=xm)
                  e=a-x;
              else
                  e=b-x;
              end;
              d=cgold*e;
          else
              d=p/q;
              u=x+d;
              if (u-a < tol2 || b-u < tol2);
                  d=tol1*sign(xm-x);
              end;
          end;
      else
          if (x>=xm)
             e=a-x;
          else
             e=b-x;
          end;
             d=cgold*e;
      end;
      if abs(d) >= tol1;
          u=x+d;
      else
          u=x+tol1*sign(d);
      end;
      fu=ls(u,xfunc,yfunc);
      if (fu <= fx);
          if (u >= x);
              a=x;
          else
              b=x;
          end;
          v=w;
          fv=fw;
          w=x;
          fw=fx;
          x=u;
          fx=fu;
      else
          if (u < x);
              a=u;
          else;
              b=u;
          end;
          if (fu <= fw || w==x);
              v=w;
              fv=fw;
              w=u;
              fw=fu;
          elseif (fu <= fv || v==x || v==w);
              v=u;
              fv=fu;
          end;
      end;
  end;
  xmin=x;
  fxmin=fx;
end