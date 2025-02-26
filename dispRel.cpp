#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <iostream>

/// Tolerance
constexpr double epsilon=
  2*std::numeric_limits<double>::epsilon();

constexpr double kappaFree=0.125;

bool twistedMass;

#define CRASH(...)							\
  internalCrash(__LINE__,__FILE__,__PRETTY_FUNCTION__,__VA_ARGS__)

using std::array;
using std::cin;
using std::cout;
using std::cerr;
using std::swap;
using std::endl;

/// Debug
__attribute__((format (printf,4,5),noreturn))
void internalCrash(const int& line,
		   const char* file,
		   const char* func,
		   const char* format,
		   ...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  cerr<<"ERROR in function "<<func<<" at line "<<line<<" of file "<<file<<": \""<<buffer<<"\""<<endl;
  exit(1);
}

/// Check if two quantities have the same sign
template <typename T>
bool sameSign(const T &a,
	      const T &b)
{
  return (a<=0 and b<=0) or (a>=0 and b>=0);
}

/// Crash if the two points do not bracket a zero
void checkNotSameSign(const double& fa,
		      const double& fb)
{
  if(sameSign(fa,fb))
    CRASH("f(a) and f(b) do not have opposite sign: %lg %lg",fa,fb);
}

/// Check if |a-b|<=|a+b|
bool areClose(const double& a,
		      const double& b,
		      const double& relTol)
{
  return fabs(a-b)<=relTol*(fabs(a+b));
}


/// Solve using Brent method
template <typename F>
double BrentSolve(F&& fun,
		  double a,
		  double b,
		  const double& relTol=epsilon*10)
{
  double fa=fun(a),fb=fun(b);
  double d=epsilon,s=0;
  
  double c=a,fc=fa;
  bool mflag=true;
  int iter=0;
  
  
  while(fabs(fa)>epsilon and
	fabs(fb)>epsilon and
	not areClose(a,b,relTol))
    {
      if(sameSign(fa,fb))
	CRASH("Something went wrong");
      
      if(not areClose(fa,fc,relTol) and
	 not areClose(fb,fc,relTol))
	{
	  s=a*fb*fc/((fa-fb)*(fa-fc))+
	    b*fc*fa/((fb-fc)*(fb-fa))+
	    c*fa*fb/((fc-fa)*(fc-fb));
	}
      else
	s=b-fb*(b-a)/(fb-fa);
      
      const double tmp=(3*a+b)/4;
      
      const double whatToComp=
	mflag?fabs(b-c):fabs(c-d);
      
      if((not ((s>tmp and s<b) or
	    (s<tmp and s>b))) or
	 fabs(s-b)>=whatToComp/2 or
	  whatToComp<epsilon)
	{
	  s=(a+b)/2;
	  mflag=1;
	}
      else
	mflag=0;
      
      const double fs=fun(s);
      d=c;
      c=b;
      fc=fb;
      
      //replace the point with opposite sign
      if(sameSign(fb,fs))
	{
	  b=s;
	  fb=fs;
	}
      else
	{
	  a=s;
	  fa=fs;
	}
      
      //makes again a the larger
      if(fabs(fa)<fabs(fb))
	{
	  swap(a,b);
	  swap(fa,fb);
	}
      iter++;
      
      if(iter>100000) CRASH("Error");
    }
  
  return b;
}

/// Compute the energy of a fermion
double tmFermionEnergy(const double& kappa,
		       const double& mu,
		       const array<double,3>& bc,
		       const double& L)
{
  const double m0=0.5/kappa-4;
  const double m2=m0*m0+mu*mu;
  
  double p2=0,p4=0;
  for(int i=0;i<3;i++)
    {
      const double p=M_PI*bc[i]/L;
      const double sinph=sin(p/2);
      const double sinph2=sinph*sinph;
      const double sinph4=sinph2*sinph2;
      p2+=sinph2;
      p4+=sinph4;
    }
  p2*=4;
  p4*=4;
  
  return 2*asinh(sqrt((m2+p2*(1+m0)+p2*p2/4-p4)/(1+m0+p2/2)/4));
}

/// Compute the energy of a naive massless fermion
double naiveMasslessFermionEnergy(const array<double,3>& bc,
				  const double& L)
{
  double sinh2E=0;
  for(int i=0;i<3;i++)
    sinh2E+=pow(sin(M_PI*bc[i]/L),2);
  
  return asinh(sqrt(sinh2E));
}

/// Returns the kappa corresponding to a given mass
double kappaOfM0(const double& m0)
{
  return 0.5/(m0+4);
}

/// Finds the mass corresponding to the external one, taking into account lattice artifacts
double mOutOfInverseDispRel(const double& M)
{
  return BrentSolve([M](const double& m)
  {
    const double kappa=
      twistedMass?
      kappaFree:
      kappaOfM0(m);
    
    const double mu=
      twistedMass?m:0.0;
    
    return tmFermionEnergy(kappa,mu,{0,0,0},1)-M;
  },
    M/2,M*2);
}

int main()
{
  cout<<"Twisted mass regularization [0-1]? ";
  cin>>twistedMass;
  
  double mMes;
  cout<<"mMes? ";
  cin>>mMes;
  
  double mLep;
  cout<<"mLep? ";
  cin>>mLep;
  
  double L;
  cout<<"L? ";
  cin>>L;
  
  int nDim;
  cout<<"NDim [1-3]? ";
  cin>>nDim;
  
  if(nDim>3 or nDim<1)
    CRASH("asked to spread momentum over %d directions, please choose in the range [1-3])",nDim);
  
  const double m=mOutOfInverseDispRel(mLep);
  
  const auto lepsEn=
    [=](const double& th)
    {
      array<double,3> bc;
      for(int i=0;i<nDim;i++)
	bc[i]=th;
      
      const double Enu=
	naiveMasslessFermionEnergy(bc,L);
      
      const double kappa=
	twistedMass?
	kappaFree:
	kappaOfM0(m);
      
      const double mu=
	twistedMass?m:0.0;
      
      const double Elep=
	tmFermionEnergy(kappa,mu,bc,L);
      
      return std::make_pair(Enu,Elep);
    };
  
  const double th=
    BrentSolve([=](const double& th)
  {
    const auto [Enu,Elep]=lepsEn(th);
    
    return mMes-Enu-Elep;
  },0.00,1000);
  
  cout<<th<<endl;
  
  const auto [Enu,Elep]=lepsEn(th);
  
  cout<<"Theta: {"<<th;
  for(int i=1;i<3;i++)
    cout<<","<<((i<nDim)?th:0.0);
  cout<<"}"<<endl;
  
  cout<<"Enu: "<<Enu<<endl;
  
  cout<<"Elep: "<<Elep<<endl;
  
  cout<<"mMes-Enu-Elep: "<<mMes-Enu-Elep<<endl;
  
  return 0;
}
