// Optional .NET Framework accelerator for india10_simulate.m.
// All historical inputs are passed by MATLAB; no fitted data are embedded.
// State layout is MATLAB column-major: eight compartments per age group.
using System;
namespace NDMC {
public static partial class India10Kernel {
    private static void Rhs(double[] x, double phase, double[] alpha,
        double[] maternal, double birth, double death, double v1, double v2,
        double sia, double[] omega, double[] maternalOmega, double waning,
        double[] dx, double[] incidence, bool sameAgeWaning) {
        double total=0, pressure=0;
        for(int i=0;i<10;i++) {
            double n=0; for(int k=0;k<8;k++) n+=x[8*i+k];
            total+=n; pressure+=x[8*i+3]/n;
        }
        double season=1-0.05*Math.Cos(2*Math.PI*phase/26);
        for(int i=0;i<10;i++) {
            int z=8*i;
            for(int k=0;k<8;k++) {
                double age=(k==0 ? maternalOmega[i] : omega[i]);
                dx[z+k]=-age*x[z+k]-death*x[z+k];
                if(i>0) dx[z+k]+=(k==0 ? maternalOmega[i-1] : omega[i-1])*x[z-8+k];
            }
            double lambda=alpha[i]*pressure*season;
            double im=maternal[i]*x[z]*lambda, si=x[z+1]*lambda, fi=x[z+2]*lambda;
            incidence[i]=im+si+fi;
            dx[z]-=im; dx[z+1]-=si; dx[z+2]-=fi;
            dx[z+3]+=incidence[i]-x[z+3]; dx[z+4]+=x[z+3];
        }
        dx[0]+=birth*total;
        if(sameAgeWaning) {
            for(int i=0;i<10;i++) { double w=waning*x[8*i]; dx[8*i]-=w; dx[8*i+1]+=w; }
        } else { double w=waning*x[16]; dx[16]-=w; dx[25]+=w; }
        double vm=v1*x[24], vs=v1*x[25];
        dx[24]-=vm; dx[25]-=vs; dx[26]+=0.05*(vm+vs); dx[29]+=0.95*(vm+vs);
        double second=0.98*v2*x[34]; dx[34]-=second; dx[38]+=second;
        for(int i=3;i<9;i++) for(int k=0;k<3;k++) {
            double v=0.70*sia*x[8*i+k]; dx[8*i+k]-=v; dx[8*i+7]+=v;
        }
    }
    public static double[] Simulate(double[] p, double[] initial, double[] mcv1,
        double[] mcv2, double[] sia, double[] omega, double[] maternalOmega,
        double waning, int substeps) {
        return SimulateCore(p,initial,mcv1,mcv2,sia,omega,maternalOmega,waning,substeps,false,0.3686);
    }
    public static double[] SimulateManuscript(double[] p, double[] initial, double[] mcv1,
        double[] mcv2, double[] sia, double[] omega, double[] maternalOmega,
        double waning, int substeps) {
        return SimulateCore(p,initial,mcv1,mcv2,sia,omega,maternalOmega,waning,substeps,true,1.0);
    }
    private static double[] SimulateCore(double[] p, double[] initial, double[] mcv1,
        double[] mcv2, double[] sia, double[] omega, double[] maternalOmega,
        double waning, int substeps, bool sameAgeWaning, double transmissionScale) {
        if(p.Length!=6 || initial.Length!=80 || mcv1.Length!=54 || mcv2.Length!=54 ||
           sia.Length!=1404 || omega.Length!=10 || maternalOmega.Length!=10 || substeps<1)
            throw new ArgumentException("Invalid India ten-age dimensions.");
        double[] x=(double[])initial.Clone(), alpha=new double[10], maternal=new double[10];
        for(int i=0;i<10;i++) alpha[i]=transmissionScale*p[i<4?0:(i<8?1:(i==8?2:3))];
        maternal[0]=0.2*p[5]; maternal[1]=0.6*p[5]; maternal[2]=p[5];
        double[] a=new double[80],b=new double[80],c=new double[80],d=new double[80];
        double[] xb=new double[80],xc=new double[80],xd=new double[80];
        double[] ia=new double[10],ib=new double[10],ic=new double[10],id=new double[10];
        double[] annual=new double[540], result=new double[24];
        double h=1.0/substeps;
        for(int step=0;step<1404;step++) {
            int year=step/26;
            double birth=(40-(24.0/52)*year)/(1000*26), death=(17-(8.0/52)*year)/(1000*26);
            for(int sub=0;sub<substeps;sub++) {
                double phase=step+1+sub*h;
                Rhs(x,phase,alpha,maternal,birth,death,mcv1[year],mcv2[year],sia[step],omega,maternalOmega,waning,a,ia,sameAgeWaning);
                for(int k=0;k<80;k++) xb[k]=x[k]+h*a[k]/2;
                Rhs(xb,phase+h/2,alpha,maternal,birth,death,mcv1[year],mcv2[year],sia[step],omega,maternalOmega,waning,b,ib,sameAgeWaning);
                for(int k=0;k<80;k++) xc[k]=x[k]+h*b[k]/2;
                Rhs(xc,phase+h/2,alpha,maternal,birth,death,mcv1[year],mcv2[year],sia[step],omega,maternalOmega,waning,c,ic,sameAgeWaning);
                for(int k=0;k<80;k++) xd[k]=x[k]+h*c[k];
                Rhs(xd,phase+h,alpha,maternal,birth,death,mcv1[year],mcv2[year],sia[step],omega,maternalOmega,waning,d,id,sameAgeWaning);
                for(int k=0;k<80;k++) {
                    x[k]+=h*(a[k]+2*b[k]+2*c[k]+d[k])/6;
                    if(Double.IsNaN(x[k]) || Double.IsInfinity(x[k]) || x[k]<-1e-10) {
                        for(int j=0;j<24;j++) result[j]=Double.NaN; return result;
                    }
                }
                for(int i=0;i<10;i++) annual[10*year+i]+=h*(ia[i]+2*ib[i]+2*ic[i]+id[i])/6;
            }
        }
        for(int j=0;j<24;j++) {
            double n=0; for(int i=0;i<10;i++) n+=annual[10*(j+30)+i];
            result[j]=p[4]*1e6*n;
        }
        return result;
    }
}
}
