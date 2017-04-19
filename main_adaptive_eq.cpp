// Finite difference method. "Equilibrium" boundary condition.
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cfloat>
#include <algorithm>    // std::min
#include <sstream>

double next_value(std::ifstream *infile)
{
  long h;
  char c, tmp[256];

  h=0;
  while(true)
  {
    c = infile->get();
    if (infile->eof()) break;
    while (c == '#')
    {
      while(c!='\n') c = infile->get();
      c = infile->get();
      if (infile->eof()) break;
    }
    if (c == '\n' || c == '\t')
    {
      tmp[h] = '\n';
      break;
    };
    tmp[h] = c;
    h++;
  };

return atof(tmp);
};

double E0AB,E0SAM,Ein,Efin,tin,v,R,F,D,T,alphaSAM,alphaAB,b,d,k0AB,k0SAM;
double ABULK,area,f,G0,t_max,t,dt,dtt;
unsigned long Nx;
//unsigned long long Nt;

using namespace std;


double E(double tt)
{
    double t1 = (Efin-Ein)/v;
    if (tt<=tin)
        return Ein;
    else if(tt<=t1+tin)
        return Ein+v*(tt-tin);
    else if(tt<=2.*t1+tin)
        return Efin+v*(t1-tt+tin);
    else
        return Ein;
};


long read_data()
{
    ifstream *ifile;
    ifile = new ifstream("data.ini");

    E0AB  = next_value(ifile);       // standard A/B potential, V
    E0SAM = next_value(ifile);

    Ein   = next_value(ifile);       // initial electrode potential (stationary), V
    Efin  = next_value(ifile);       // largest potential value, V

    tin   = next_value(ifile);
    v     = next_value(ifile);        // potential cycle speed, V/s
    R     = next_value(ifile);  // gas constant, J/(K*mol)
    F     = next_value(ifile); // Faraday constant, C
    D     = next_value(ifile);       // diffusion coefficient, cm2/sec
    T     = next_value(ifile);   // temperature, K
    alphaSAM = next_value(ifile);        // alpha for SAM
    alphaAB= next_value(ifile);        // alpha for A/B
    b     = next_value(ifile);          // beta, dimensionless
    d     = next_value(ifile);          // distance from el to SAM in Fc radiuses (dimensionless)

    k0AB  = next_value(ifile);      // rate constant for A/B, cm/sec
    Nx   = next_value(ifile);
    dtt   = next_value(ifile);

    k0SAM = next_value(ifile);        // rate constant for SAM, 1/sec

    f     = F/(R*T);

    ABULK = next_value(ifile);     // bulk concentration of A, mol/cm3
    area  = next_value(ifile);     // surface area, cm2
    G0  = next_value(ifile);      // surface conc of SAM, mol/cm2
    t_max    = tin+2.*(Efin-Ein)/v;  //

    ifile->close();
    delete ifile;

    return 0;
};


long push_inidata(std::stringstream &ofile)
{
    ofile << "E0AB = " << E0AB << endl <<
    "E0SAM = " << E0SAM << endl <<
    "Ein = " << Ein << endl <<
    "Efin = " << Efin << endl <<
    "tin = " << tin << endl<<
    "v = " << v << endl <<
    "R = " << R << endl <<
    "F = " << F << endl <<
    "D = " << D << endl <<
    "T = " << T << endl <<
    "alphaSAM = " << alphaSAM << endl <<
    "alphaAB = " << alphaAB << endl <<
    "b = " << b << endl <<
    "d = " << d << endl <<
    "k0AB = " << k0AB << endl <<
    "dt accuracy = " << dtt << endl <<
    "k0SAM = " << k0SAM << endl <<
    "ABULK = " << ABULK << endl <<
    "area = " << area << endl <<
    "Gamma0 = " << G0 << endl <<
    "t = " << t_max << endl;

    return 0;
}


time_t print_time_stamp()
{
    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    cout<<asctime(timeinfo);

    return rawtime;
};


////////////////////////////////////////////////////////////////////////
int main( int argc, char *argv[] ) {

    double *A, *AP, *Atmp, Gox=0, GoxP=0, dx=0;
    double x_max=0, xi=0, kSAM=0, k=0;

    time_t time_start, time_finish;

    read_data();

    x_max = 6.*sqrt(D*t_max);
    dx = x_max/(double)Nx;

    cout<<"Time grid: "<<endl<<"-------------------------------"<<endl;
    cout<<"dt, s    = "<<dtt<<endl;
    cout<<"t_max, s = "<<t_max<<endl;
    //cout<<"Nt       = "<<Nt<<endl<<endl;
    cout<<"X grid:     "<<endl<<"-------------------------------"<<endl;
    cout<<"dx, cm    = "<<dx<<endl;
    cout<<"x_max, cm = "<<x_max<<endl;
    cout<<"Nx        = "<<Nx<<endl<<endl;

    //ofstream *ofile = new ofstream("out.dat");
    //ofstream *Cfile = new ofstream("C.dat");

    std::stringstream ofile, Cfile;

    push_inidata(ofile);

    ofile<<"dt\ttime\tE(t)\tEsam(t)\tcurr_e+s\tcurr_el"
    <<"\t"<<"[A]_0"<<endl;

    // allocate concentration
    A  = new double[Nx+1];
    AP = new double[Nx+1];

    // initial conditions
    for (unsigned long j=0;j<=Nx; j++) A[j] = ABULK;
    Gox = G0/(1.+exp(-f*(E(0)-E0SAM)));

    dt = 1e-4;
    t = 0;

    // iterating in time
    cout<<"Starting at: ";
    time_start = print_time_stamp();
    cout<<endl;


    double I_el=0, I_sam=0;


    //main loop

    cout << scientific << setprecision(2);

    for(double s=1e-3; s<=1; s+=1e-3)
    {
        // saving currents etc.
        ofile << dt<<"\t"<<t<<"\t"<<E(t)<<"\t"
                <<E0SAM+log(Gox/(G0-Gox))/f
                <<"\t"<<I_sam*area*F
                <<"\t"<<I_el*area*F
                <<"\t"<<A[0]<<endl;

        // saving concentration profiles
        for (unsigned long j=0; j<Nx; j++)
        {
            Cfile << A[j] << "\t";
        };
        Cfile << A[Nx] << endl;


        time_t rawtime;
        int hrs, min;

        time (&rawtime);
        rawtime -= time_start;
        hrs = rawtime/3600;
        rawtime %= 3600;
        min = rawtime/60;
        rawtime %= 60;

        cout << setfill(' ') << setw(2) << (long)(t/t_max*100) << "%" << "    ";
        cout << "E=" << E(t)
             << " I=" << I_el * area * F
             << " IS=" << I_sam * area * F
             << " kS=" << kSAM
//             << " Gox=" << Gox
             << " dt=" << dt
             << " (" << setfill('0') << setw(2) << hrs << ":"
             << setw(2) << min << ":"
             << setw(2) << rawtime << ")" << endl;

        while(t<s*t_max)
        {
            // 1) find Gox+
            xi = f*(E(t) - E0SAM);
            kSAM = k0SAM*exp(-alphaSAM*xi)*exp(-b*d);
            I_el = kSAM*exp(xi)*G0-kSAM*(1+exp(xi))*Gox;
            I_sam = I_el-D*(A[1]-A[0])/dx;

            // evaluate dt

            //tmp = dtt/(k0SAM*exp(alphaSAM*abs(xi))*exp(-b*d));
            //tmp = 3e-20/abs(I_sam);

            //find dt from eqn for Gox
            double dtG = abs(Gox/I_sam) * dtt;
            if(dtG > 1e-6) dtG = 1e-6;

            //find dt from Aj
            double dtA = 1e0;
            for (unsigned long j=1; j<Nx; j++)
            {
                dtA = std::min(dtA, abs(A[j]/(A[j-1] - 2.*A[j] + A[j+1])) *
                               (pow(dx,2)/D) * 1e-10);
            };


            dt = std::min(dtA, dtG);

            //if(dt==dtA) std::cout<<dtA<<std::endl;

            GoxP = Gox + dt*I_sam;

            // 2) find Aj+
            for (unsigned long j=1; j<Nx; j++)
            {
                AP[j] = A[j] + D*dt/pow(dx,2)*(A[j-1] - 2.*A[j] + A[j+1]);
            };
            AP[Nx] = ABULK;


            // 3) find A0+
            k = exp(f*(E0AB-E0SAM));
            AP[0] = k*(G0-Gox)/(Gox+k*(G0-Gox))*ABULK;

            // update values
            Gox = GoxP;

            Atmp = A;
            A = AP;
            AP = Atmp;

            t += dt;
        };
    };

    cout<<"Finished at: ";
    time_finish = print_time_stamp();
    cout<<endl<<endl;
    cout<<"Needed "<<difftime(time_finish, time_start)<<" s";

    cout<<endl;

    // flushing data to disk

    ofstream outfile("out.dat");
    outfile << ofile.str();
    outfile.close();

    ofstream concfile("C.dat");
    concfile << Cfile.str();
    concfile.close();

    //ofile->close();
    //Cfile->close();

    //delete ofile;
    //delete Cfile;
    delete [] A ;
    delete [] AP;

    return 0;
};

