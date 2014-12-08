
#include "DetectorModel.h"

#ifdef HAVE_ROOT
#include "Rtypes.h"
// #include "TFile.h"
#endif

#include <cstdio>
#include <iostream>

// #define DEBUG

using namespace DetectorModel;

//////////////////////
// Static constants //
//////////////////////

// -- units
double Units::GeV = 1.0;
double Units::MeV = 1./1000.;
double Units::cm  = 10.;
double Units::mm  = 1.;

// -- hadronic grid
size_t Defaults::m_nEtaBinsHAD = 100;
double Defaults::m_etaMinHAD   = -5.;
double Defaults::m_etaMaxHAD   =  5.;
size_t Defaults::m_nPhiBinsHAD =  64; 

// -- electromagnetic grid
size_t Defaults::m_nEtaBinsEM  = 100;
double Defaults::m_etaMinEM    = -2.5;
double Defaults::m_etaMaxEM    =  2.5;
size_t Defaults::m_nPhiBinsEM  = 256; 

// -- calorimeter dimensions
double Defaults::m_caloRadius    = 1200.*Units::mm;
double Defaults::m_caloHalfWidth = 2500.*Units::mm;
double Defaults::m_caloEtaMin    = -5.;
double Defaults::m_caloEtaMax    =  5.;

///////////////////////////
// Useful tool functions //
///////////////////////////

// -- distance between two points
double Tools::distance(const SpacePoint& pointA,const SpacePoint& pointB)
{
  double dx(pointA.x()-pointB.x());
  double dy(pointA.y()-pointB.y());
  double dz(pointA.z()-pointB.z());
  return sqrt(dx*dx+dy*dy+dz*dz);
}

// -- (shortest) distance between a line and a point 
double Tools::distance(const Line& line,const SpacePoint& point)
{
  double t(parallel(line,point));
  double dx(point.x()-line.origin().x()-line.dirX()*t);
  double dy(point.y()-line.origin().y()-line.dirY()*t);
  double dz(point.z()-line.origin().z()-line.dirZ()*t);
  return std::sqrt(dx*dx+dy*dy+dz*dz);
}

//////////////////
//
// --------------------------------------------- Geometry: Point
//

SpacePoint::SpacePoint() : m_x(0.), m_y(0.), m_z(0.), m_rad(0.)
{ }

SpacePoint::SpacePoint(double x,double y,double z)
  : m_x(x), m_y(y), m_z(z)
{
  m_rad = std::sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
}

SpacePoint::SpacePoint(const SpacePoint& point)
  : m_x(point.x()), m_y(point.y()), m_z(point.z())
{
  m_rad = std::sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
}

SpacePoint::~SpacePoint()
{ }

//
// --------------------------------------------- Geometry: Line
//

Line::Line() 
  : m_origin(0.,0.,0.), m_dirX(0.), m_dirY(0.), m_dirZ(1.0), m_length(0.)
{ }

Line::Line(const Line& line)
  : m_origin(line.origin()), m_dirX(line.dirX()), m_dirY(line.dirY()),
    m_dirZ(line.dirZ()), m_length(line.length())
{ }

Line::Line(const SpacePoint& a,const SpacePoint& b)
  : m_origin(a)
{
  double r(a.distance(b));
  if ( r > 0 )
    {
      m_dirX   = (b.x()-m_origin.x())/r;
      m_dirY   = (b.y()-m_origin.y())/r;
      m_dirZ   = (b.z()-m_origin.z())/r;
      m_length = r;
    }
  else
    {
      m_dirX   = 0.;
      m_dirY   = 0.;
      m_dirZ   = 1.;
      m_length = 0.;
    }
} 

Line::~Line()
{ }

//
// --------------------------------------------- Profiles: tower grid
//

Grid::Grid() 
  : m_grid(0)
  , m_nBins(0)
  , m_nEta(0)
  , m_nPhi(0) 
  , m_etaMin(0.)
  , m_etaMax(0.)
  , m_dEta(0.)
  , m_dPhi(0.)
  , m_twoPi(4.*std::asin(1.))
{ }

Grid::Grid(size_t nEta,double etaMin,double etaMax,size_t nPhi)
  : m_grid(nEta*nPhi,0.)
  , m_nBins(nEta*nPhi)
  , m_nEta(nEta)
  , m_nPhi(nPhi)
  , m_etaMin(etaMin)
  , m_etaMax(etaMax)
  , m_dEta(0.)
  , m_dPhi(0.)
  , m_twoPi(4.*std::asin(1.))
{
  m_dEta = m_nEta != 0 
    ? ( m_etaMax - m_etaMin ) / static_cast<double>(m_nEta)
    : 0.;
  m_dPhi = m_nPhi != 0
    ?  m_twoPi / static_cast<double>(m_nPhi)
    : 0.;
}

Grid::Grid(const Grid& grid)
  : m_grid(grid.m_grid)
  , m_nBins(grid.m_nBins)
  , m_nEta(grid.m_nEta)
  , m_nPhi(grid.m_nPhi)
  , m_etaMin(grid.m_etaMin)
  , m_etaMax(grid.m_etaMax)
  , m_dEta(grid.m_dEta)
  , m_dPhi(grid.m_dPhi)
  , m_twoPi(grid.m_twoPi)
{ }

Grid::~Grid()
{ }

bool Grid::etaIndex(double eta,index_t& index) const
{
  if ( eta >= m_etaMin && eta < m_etaMax ) 
    {
      index = static_cast<index_t>(std::floor((eta-m_etaMin)/m_dEta));
      return true;
    }
  return false;
} 

bool Grid::phiIndex(double phi,index_t& index) const
{
  if ( phi < 0. ) phi += m_twoPi;
  index = static_cast<index_t>(std::floor(phi/m_dPhi));
  while ( index >= m_nPhi ) { index -= m_nPhi; }
  return true;
}

particle_t Grid::pseudoJet(index_t index) const
{
  if ( index < m_nBins )
    {
      double e(this->energy(index));
      double theta(std::acos(std::tanh(this->eta(index))));
      double phi(this->phi(index));
      return particle_t(e*std::sin(theta)*std::cos(phi),
			e*std::sin(theta)*std::sin(phi),
			e*std::cos(theta),e);
    }
  else
    {
      return particle_t(0.,0.,0.,0.);
    }
}

container_t Grid::fillNoZero() const
{
  std::vector<double>::const_iterator fGrid(m_grid.begin());
  std::vector<double>::const_iterator lGrid(m_grid.end());
  container_t jetList(0);
  index_t index(0);
  for ( ; fGrid != lGrid; ++fGrid, index++ )
    {
      if ( *fGrid > 0. ) 
	jetList.push_back(this->pseudoJet(*fGrid,
					  this->eta(index),
					  this->phi(index)));
    }
  return jetList;
}

container_t Grid::fillAll() const
{
  std::vector<double>::const_iterator fGrid(m_grid.begin());
  std::vector<double>::const_iterator lGrid(m_grid.end());
  container_t jetList(0);
  index_t index(0);
  for ( ; fGrid != lGrid; ++fGrid, index++ )
    {
      jetList.push_back(this->pseudoJet(*fGrid,
					this->eta(index),
					this->phi(index)));
    }
  return jetList;
}

bool Grid::add(double px,double py,double pz,double e)
{
  // check ranges;
  double pt(px*px+py*py);
  if ( pt == 0. ) return false;
  if ( pz == 0. ) pz = pt/1000.;
  double p(sqrt(pt+pz*pz));

  pt = sqrt(pt);
  double eta(-0.5*log((p+pz)/(p-pz)));
  double phi(std::atan2(py,px));
  return this->add(e,eta,phi);
} 

//
// --------------------------------------------- Detector description
//

DetectorDescription::DetectorDescription()
{ }

DetectorDescription::~DetectorDescription()
{ }

bool DetectorDescription::inAcceptance(const particle_t& particle) const
{ 
  return this->inAcceptance(particle.e(),particle.pseudorapidity(),
			    particle.phi(),particle.m());
}


SpacePoint DetectorDescription::impact(const particle_t& particle) const
{ 
  return this->impact(particle.e(),particle.pseudorapidity(),particle.phi(),
		      particle.m()); 
}

CylindricalCalo::CylindricalCalo(double r,double z,double etaMin,double etaMax)
  : DetectorDescription()
  , m_radius(r)
  , m_width(z)
  , m_etaMin(etaMin)
  , m_etaMax(etaMax)
{
  double theta(std::atan(r/z));
  double eta(-std::log(std::tan(theta/2.)));
  m_etaBoundFwd = fabs(eta);
  m_etaBoundBwd = -fabs(eta);
 #ifdef DEBUG
  std::cout << "[CylindricalCalo] constructed with:" << std::endl;
  std::cout << "                  r      = " << m_radius << std::endl;
  std::cout << "                  dz     = " << m_width  << std::endl;
  std::cout << "                  theta  = " << theta << std::endl;
  std::cout << "                  etaFwd = " << m_etaBoundFwd << std::endl;
  std::cout << "                  etaBck = " << m_etaBoundBwd << std::endl;
  std::cout << "                  etaMin = " << m_etaMin << std::endl;
  std::cout << "                  etaMax = " << m_etaMax << std::endl;
 #endif
}

CylindricalCalo::~CylindricalCalo()
{ }

DetectorDescription* CylindricalCalo::clone() const
{
  return new CylindricalCalo(m_radius,m_width,m_etaMin,m_etaMax);
}

SpacePoint CylindricalCalo::impact(double e,double eta,double phi,double m)
 const
{
  // where are we 
  if ( !this->inAcceptance(e,eta,phi,m) ) return SpacePoint();

  // radial part
  double theta(std::acos(std::tanh(eta)));
  if ( eta > m_etaBoundBwd && eta < m_etaBoundFwd )
    {
      double rVtx(m_radius/std::sin(theta));
      double x(m_radius*std::cos(phi));
      double y(m_radius*std::sin(phi));
      double z(rVtx*std::cos(theta));
      return SpacePoint(x,y,z);
    }
  // longitudinal part
  else
    {
      double rPerp(m_width*std::tan(theta));
      double x(rPerp*std::cos(phi));
      double y(rPerp*std::sin(phi));
      double z = eta < 0. ? -m_width : m_width;
      return SpacePoint(x,y,z);
    }
}
//
// --------------------------------------------- Shower shapes
//
#ifndef HAVE_ROOT
RadialShape::RadialShape(double width) : m_extension(width) { }
#else
RadialShape::RadialShape(double width,const std::string& name) 
  : m_extension(width) 
  , m_name(name)
  , m_histBooked(false)
{ }
#endif

RadialShape::~RadialShape() { }

double RadialShape::evaluate(double r,double e)
{
  return (this->operator())(r,e);
}

double RadialGaussShape::m_probLvl = 4.;

#ifndef HAVE_ROOT
RadialGaussShape::RadialGaussShape(double width) 
  : RadialShape(m_probLvl*width)
#else
RadialGaussShape::RadialGaussShape(double width,const std::string& name) 
    : RadialShape(m_probLvl*width,name)
    , d_profile((TH2D*)0)
#endif
{
  m_width = width;
  //  this->setExtension(m_probLvl*m_width);
  this->setParameters();
  std::cout << "[RadialGaussShape] sigma/extend: " << m_width << "/" << this->extension() << std::endl;
}
#ifndef HAVE_ROOT
RadialGaussShape::RadialGaussShape(double width,double rMax)
  : RadialShape(rMax)
#else
RadialGaussShape::RadialGaussShape(double width,double rMax,const std::string& name)
  : RadialShape(rMax)
  , d_profile((TH2D*)0)
#endif
{ 
  m_width = width;
  // this->setExtension(rMax);
  this->setParameters();
}

RadialGaussShape::RadialGaussShape(const RadialGaussShape& shape)
#ifndef HAVE_ROOT
  : RadialShape(shape.extension())
#else
  : RadialShape(shape.extension(),shape.name())
  , d_profile((TH2D*)0)
#endif
  , m_width(shape.m_width)
  , m_const(shape.m_const)
  , m_width2(shape.m_width2)
{
  this->setParameters();
}

RadialGaussShape::~RadialGaussShape()
{ }

void RadialGaussShape::setParameters()
{
  m_width2 = 2.*m_width*m_width;
  m_const  = 1.; // / std::sqrt(fabs(4*std::asin(1.))*m_width2);
  this->normalize();
}

void RadialGaussShape::normalize()
{
  double r(0.);
  double dr(this->extension()/100.);
  double norm(0.);
  while ( r < this->extension() )
    {
      double f((this->operator())(r));
      norm += f;
      //      std::cout << " .... r/f " << r << "/" << f << std::endl;
      r += dr;
    }
  m_const = 0.5/norm;
  std::cout << "[RadialGaussShape::normalize()] norm/const " << norm << "/" << m_const << std::endl;
}

double RadialGaussShape::operator()(double r,double /*e*/) const
{
#ifndef HAVE_ROOT
  if ( r > this->extension() ) return 0.;
  return m_const * std::exp(-r*r/m_width2);
#else
  double f = r > this->extension() ? 0. : m_const * std::exp(-r*r/m_width2);
  if ( m_histBooked ) d_profile->Fill(r,std::log10(f));
  return f;
#endif
}

void RadialGaussShape::print() const
{
  std::cout << "[RadialGaussShape] (constant,width) = (" << m_const << "," << m_width << ")" << std::endl;
  double r(-this->extension());
  double dr(fabs(r)/100.);
  for ( ; r < this->extension(); r+=dr )
    {
      std::cout << " ----- f(" << r << ")*1000 = " << (this->operator())(r)*1000. << std::endl; 
    }
}

RadialShape* RadialGaussShape::clone() const
{
  return new RadialGaussShape(*this);
}
#ifdef HAVE_ROOT
bool RadialGaussShape::setupHists()
{
  //
  if ( !m_histBooked )
    {
      std::string fullName = "Profiles_" + this->name();
      double xmin(0.);
      double xmax(1.1*this->extension());
      double ymin(-4.25);
      double ymax(-0.25);
      int    nbins(100);
      d_profile = new TH2D(fullName.c_str(),fullName.c_str(),nbins,xmin,xmax,nbins,ymin,ymax);
      m_histBooked = true;
      return true;
    }
  return false;
}

bool RadialGaussShape::writeHists()
{
  if ( m_histBooked )
    {
      printf("[RadialGaussShape::writeHists()] - write histogram \042%s\042 to file\n",d_profile->GetTitle());
      //      std::string fName = "RadialGaussShape_" + this->name() + ".root";
      //      TFile* f = new TFile(fName.c_str(),"RECREATE");
      d_profile->Write();
      //      f->Close();
    }
  return true;
}
#endif

////////////////////////////////
// HAD Shower parametrization //
////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
// Hadronic shower profiles exhibit more energy dependence. Two profiles for 30 and 120
// GeV pions are available from H1. The difference between these profiles is used
// to extrapolate from low to high energy in log(E). E<30 GeV (> 210 GeV) uses the
// low (high) energy profile.
// 

// energy extrapolation 
double RadialExpShapeHAD::m_e0    = 30.*Units::GeV;
double RadialExpShapeHAD::m_e1    = 120.*Units::GeV;
double RadialExpShapeHAD::m_slope = 0.928;

// valid energy range
double RadialExpShapeHAD::m_shapeLimitLo = m_e0;
double RadialExpShapeHAD::m_shapeLimitHi = 200.*Units::GeV;

// parameterization of shape difference -------------------------------
double RadialExpShapeHAD::m_r0    = 101.*Units::mm; 
double RadialExpShapeHAD::m_r1    = 201.*Units::mm;
// ------------------------------------------------------ low range
double RadialExpShapeHAD::m_p00   =  1.03614;
double RadialExpShapeHAD::m_p01   = -0.00330184;
double RadialExpShapeHAD::m_p02   =  4.88754e-06;
double RadialExpShapeHAD::m_p03   =  4.88734e-07;
// ------------------------------------------------------ mid range
double RadialExpShapeHAD::m_p10   = -2.17745;
double RadialExpShapeHAD::m_p11   =  0.0592228;
double RadialExpShapeHAD::m_p12   = -0.000300749;
double RadialExpShapeHAD::m_p13   =  4.84098e-07;
// ------------------------------------------------------ high range
double RadialExpShapeHAD::m_p20   =  2.19605;
double RadialExpShapeHAD::m_p21   = -0.00383414;
double RadialExpShapeHAD::m_p22   =  2.12902e-06;
// parameterization of reference shape --------------------------------
double RadialExpShapeHAD::m_refp0 =  6.19709e-03;
double RadialExpShapeHAD::m_refp1 = -1.80291e-02;
double RadialExpShapeHAD::m_refp2 =  3.46430e-01;
double RadialExpShapeHAD::m_refp3 = -6.61100e-02;
// parametrization of range -------------------------------------------
double RadialExpShapeHAD::m_rangep0 =  113.48;
double RadialExpShapeHAD::m_rangep1 = -0.306114;
double RadialExpShapeHAD::m_rangep2 =  0.00138404;
double RadialExpShapeHAD::m_rangep3 = -2.39028e-06;
double RadialExpShapeHAD::m_minRange = 80.;

#ifndef HAVE_ROOT
RadialExpShapeHAD::RadialExpShapeHAD(double width) 
  : RadialShape(width)
#else
RadialExpShapeHAD::RadialExpShapeHAD(double width,const std::string& name) 
  : RadialShape(width,name)
  , d_profile((TH2D*)0)
  , d_range((TH2D*)0)
#endif
  , m_logE0(std::log(m_e0))
  , m_logE1(std::log(m_e1))
{
  this->setParameters();
}

RadialExpShapeHAD::RadialExpShapeHAD(const RadialExpShapeHAD& shape)
#ifndef HAVE_ROOT
  : RadialShape(shape.extension())
#else
  : RadialShape(shape.extension(),shape.name())
  , d_profile((TH2D*)0)
  , d_range((TH2D*)0)
#endif
  , m_logE0(shape.m_logE0)
  , m_logE1(shape.m_logE1)
{
  this->setParameters();
}

RadialExpShapeHAD::~RadialExpShapeHAD()
{ }

void RadialExpShapeHAD::setParameters()
{
  m_kappa = m_logE1/m_logE0;
  m_appak = 1. - m_kappa;
  // functional parameters for shape difference 
  m_diff1.clear();
  m_diff1.push_back(m_p00);
  m_diff1.push_back(m_p01);
  m_diff1.push_back(m_p02);
  m_diff1.push_back(m_p03);
  m_diff2.clear();
  m_diff2.push_back(m_p10);
  m_diff2.push_back(m_p11);
  m_diff2.push_back(m_p12);
  m_diff2.push_back(m_p13);
  m_diff3.clear();
  m_diff3.push_back(m_p20);
  m_diff3.push_back(m_p21);
  m_diff3.push_back(m_p22);  
  // functional parameters for range/width
  m_range.clear();
  m_range.push_back(m_rangep0);
  m_range.push_back(m_rangep1);
  m_range.push_back(m_rangep2);
  m_range.push_back(m_rangep3);
}

RadialShape* RadialExpShapeHAD::clone() const 
{
  return new RadialExpShapeHAD(*this);
}

void RadialExpShapeHAD::print() const
{
  std::cout << "[RadialExpShapeHAD] - simple parameterization of hadronic"
	    << " shower shapes in the energy rsange " 
	    << m_e0/Units::GeV << " - "
	    << m_e1/Units::GeV << " GeV (fixed outside of this range)" << std::endl;
  printf("[RadialExpShapeHAD] - principal parameters:\n");
  printf("     Rmax = %7.3f mm\n",this->extension());
  printf("     Approximate containment radii for E = (%5.1f,%5.1f) GeV\n",
	 m_e0/Units::GeV,m_e1/Units::GeV);
  printf("     r(intE/E = 50.0%) = (%7.3f,%7.3f) mm\n",
	 this->containment(0.50,m_e0),this->containment(0.50,m_e1));
  printf("     r(intE/E = 75.0%) = (%7.3f,%7.3f) mm\n",
	 this->containment(0.75,m_e0),this->containment(0.75,m_e1));
  printf("     r(intE/E = 90.0%) = (%7.3f,%7.3f) mm\n",
	 this->containment(0.90,m_e0),this->containment(0.90,m_e1));
  printf("     r(intE/E = 98.0%) = (%7.3f,%7.3f) mm\n",
	 this->containment(0.98,m_e0),this->containment(0.98,m_e1));
  printf("     r(intE/E = 99.5%) = (%7.3f,%7.3f) mm\n",
	 this->containment(0.995,m_e0),this->containment(0.995,m_e1));
}

double RadialExpShapeHAD::containment(double f,double e) const
{
  // get total integral
  double rMax(this->rangeLimit(e));
  std::cout << "Rangelimit " << rMax << " (" << e << " GeV" << std::endl;
  double ir(this->integral(rMax,e));
  std::cout << "Integral   " << ir << std::endl;
  double dr(rMax/1000.);
  double r(0.);
  double s(0.);
  while ( s/ir < f && r < rMax ) { s += (this->operator())(r,e); r += dr; }
  return r;
}

double RadialExpShapeHAD::integral(double r,double e) const
{
  double rMax(std::min(r,this->rangeLimit(e)));
  double dr(this->extension()/1000.);
  double ir(0.);
  for (double r(0.); r < rMax; r+=dr ) { ir += (this->operator())(r,e); }
  return ir;
}

double RadialExpShapeHAD::evaluate(double r,double e) const
{
  double rl(this->rangeLimit(e));
  //  double ir(this->integral(rl,e));
  double ff = r > rl ? 0. : this->refShape(r)*this->eParm(e,r);///ir;
#ifdef HAVE_ROOT
  if ( m_histBooked )
    {
      if ( ff > 0. ) d_profile->Fill(r,std::log10(ff));
      if ( e > 0.  ) d_range->Fill(std::log10(e),rl);
    }
#endif
  return ff;
}

double RadialExpShapeHAD::operator()(double r,double e) const 
{
  return this->evaluate(r,e);
}

double RadialExpShapeHAD::rangeLimit(double e) const
{
  return std::max(Math::poly(m_range,e),m_minRange);
}

double RadialExpShapeHAD::refShape(double r) const
{
  return m_refp0*std::exp(m_refp1*r) + m_refp2*std::exp(m_refp3*r);
}

double RadialExpShapeHAD::sDiff(double r) const
{
  if ( r < m_r0 ) 
    {
      return Math::poly(m_diff1,r);
    } 
  else if ( r >= m_r0 && r < m_r1 )
    {
      return Math::poly(m_diff2,r);
    }
  else
    {
      return Math::poly(m_diff3,r);
    }
}

//////////////////
// Simple log(E) interpolation between two measured shapes. The new shape
// is given given by refShape(..)/eParm(...).
//////////////////

double RadialExpShapeHAD::eParm(double e,double r) const
{
  double energy = e < m_e0 ? m_shapeLimitLo : e > m_e1 ? m_shapeLimitHi : e;
  double a((m_slope*this->sDiff(r)-m_kappa)/m_appak);
  double b((1-a)/m_logE0);
  return a + b * std::log10(energy);
}

#ifdef HAVE_ROOT
bool RadialExpShapeHAD::setupHists()
{
  if ( !m_histBooked )
    {
      printf("[RadialExpShapeHAD::setupHists(%s)] - book histogram ",this->name().c_str());
      std::string fullName = "Profiles_Energy_" + this->name();
      printf("\042%s\042\n",fullName.c_str());
      int nbins(100);
      double xmin(0.);
      double xmax(1.1*std::max(this->rangeLimit(m_e0),this->rangeLimit(m_e1)));
      double ymin(-4.25);
      double ymax(-0.25);
      d_profile = new TH2D(fullName.c_str(),fullName.c_str(),nbins,xmin,xmax,nbins,ymin,ymax);
      printf("[RadialExpShapeHAD::setupHists(%s)] - book histogram ",this->name().c_str());
      std::string rngeName = "Profiles_Range_" + this->name();
      printf("\042%s\042\n",rngeName.c_str());
      double emin(10.*Units::MeV);
      double emax(1000.*Units::GeV);
      double rmin(0.9*this->rangeLimit(emax));
      double rmax(1.1*this->rangeLimit(emin));
      d_range      = new TH2D(rngeName.c_str(),rngeName.c_str(),50,-1.,4.,nbins,rmin,rmax);
      m_histBooked = true;
      return true;
    }
  return false;
}

bool RadialExpShapeHAD::writeHists()
{
  if ( m_histBooked )
    {
      printf("[RadialExpShapeHAD::writeHists()] - write histogram \042%s\042 to file\n",d_profile->GetTitle());
      //     std::string fName = "RadiaLExpShapeHAD_" + this->name()+".root";
      //     TFile* f = new TFile(fName.c_str(),"RECREATE");
      d_profile->Write();
      printf("[RadialExpShapeHAD::writeHists()] - write histogram \042%s\042 to file\n",d_range->GetTitle());
      d_range->Write();
      //     f->Close();
    }
  return m_histBooked;
}
#endif

///////////////////////////////
// EM Shower parametrization //
///////////////////////////////

// dE/E/dr = p0 * exp(p1*r) + p2 * exp(p3*r)

double RadialExpShapeEM::m_p0 =  6.65736e-02; 
double RadialExpShapeEM::m_p1 = -1.11378e-01;///Units::cm; 
double RadialExpShapeEM::m_p2 =  3.50344e-04;
double RadialExpShapeEM::m_p3 = -2.77749e-02;///Units::cm;

#ifndef HAVE_ROOT
RadialExpShapeEM::RadialExpShapeEM(double width) 
  : RadialShape(width)
#else
  RadialExpShapeEM::RadialExpShapeEM(double width,const std::string& name) 
    : RadialShape(width,name)
    , d_profile((TH2D*)0)
#endif
  , m_norm(1.)
{
  this->normalize();
}
RadialExpShapeEM::RadialExpShapeEM(const RadialExpShapeEM& shape)
#ifndef HAVE_ROOT
  : RadialShape(shape.extension())
#else
  : RadialShape(shape.extension(),shape.name())
  , d_profile((TH2D*)0)
#endif
  , m_norm(shape.m_norm)
{ }

RadialExpShapeEM::~RadialExpShapeEM()
{ }

void RadialExpShapeEM::normalize() 
{
  m_norm = 1.0 / this->integral(this->extension());
}

double RadialExpShapeEM::operator()(double r,double /*e*/) const
{
#ifndef HAVE_ROOT 
  return this->eval(r) * m_norm;
#else
  double f(this->eval(r)*m_norm);
  if ( m_histBooked && f > 0. )
    {
      double x(r);
      double y(std::log10(f));
      d_profile->Fill(x,y);
    }
  return f;
#endif 
}

RadialShape* RadialExpShapeEM::clone() const
{
  return new RadialExpShapeEM(*this);
}

void RadialExpShapeEM::print() const
{
  std::cout << "[RadialExpShapeEM] - simple parameterization of an electromagnetic"
	    << " shower shape at fixed energy (from H1, E = 30 GeV)." << std::endl;
  std::cout << "[RadialExpShapeEM] - principal parameters: " << std::endl;
  std::cout << "                     1/E dE/(rdr) = "
	    << 1./m_norm << " * ["
	    << m_p0 << " * exp(" << m_p1 << " * r) + "
	    << m_p2 << " * exp(" << m_p3 << " * r) ]" << std::endl;
  std::cout << "                     Rmax = " << this->extension() << " mm" << std::endl;
  std::cout << "                     Approximate containment radii:" << std::endl;
  std::cout << "                     r(dE/E = 50.0%) = " << this->containment(0.5)
	    << " mm" << std::endl;
  std::cout << "                     r(intE/E = 75.0%) = " << this->containment(0.75)
	    << " mm" << std::endl;
  std::cout << "                     r(intE/E = 90.0%) = " << this->containment(0.9)
	    << " mm" << std::endl;
  std::cout << "                     r(intE/E = 98.0%) = " << this->containment(0.98)
	    << " mm" << std::endl;
  std::cout << "                     r(intE/E = 99.5%) = " << this->containment(0.995)
	    << " mm" << std::endl;
}

double RadialExpShapeEM::containment(double fraction) const
{
  double rMax(this->extension());
  double dr(rMax/100.);
  double r(0.);
  double af(fraction/m_norm);
  double integral(0.);
  while ( r < rMax && integral < af ) { integral += this->eval(r); r += dr; }
  return r;
}

double RadialExpShapeEM::exp0(double r) const
{
  double f = m_p0*std::exp(m_p1*std::abs(r));
 #ifdef DEBUG
  printf("[RadialExpShapeEM::exp0(%f)] f(%f) = %f exp(%f %f) = %f\n",r,r,m_p0,m_p1,r,f);
#endif
  return f;
}

double RadialExpShapeEM::exp1(double r) const
{
  double f = m_p2*std::exp(m_p3*r);
 #ifdef DEBUG
  printf("[RadialExpShapeEM::exp1(%f)] f(%f) = %f exp(%f %f) = %f\n",r,r,m_p2,m_p3,r,f);;
#endif
  return f;
}

#ifdef HAVE_ROOT
bool RadialExpShapeEM::setupHists()
{
  if ( !m_histBooked )
    {
      printf("[RadialExpShapeEM::setupHists(%s)] - book histogram ",this->name().c_str());
      std::string fullName = "Profiles_" + this->name();
      printf("\042%s\042\n",fullName.c_str());
      int nbins(100);
      double xmin(0.);
      double xmax(1.1*this->extension());
      double ymin(-4.25);
      double ymax(-0.25);
      d_profile = new TH2D(fullName.c_str(),fullName.c_str(),nbins,xmin,xmax,nbins,ymin,ymax);
      m_histBooked = true;
      return true;
    }
  return false;
}

bool RadialExpShapeEM::writeHists()
{
  if ( m_histBooked )
    {
      printf("[RadialExpShapeEM::writeHists()] - write histogram \042%s\042 to file\n",d_profile->GetTitle());
      //      std::string fName = "RadialExpShapeEM_" + this->name() + ".root";
      //      TFile* f = new TFile(fName.c_str(),"RECREATE");
      d_profile->Write();
      //     f->Close();
    }
  return m_histBooked;
}
#endif

double RadialProfile::m_relPrec = 1./1000.;

RadialProfile::RadialProfile(double rPerp,
			     size_t nWidthDiv, 
			     size_t nPhiDiv,
			     const RadialShape& shape)
  : m_shape(shape.clone())
  , m_widthDiv(nWidthDiv)
  , m_phiDiv(nPhiDiv)
{
#ifdef HAVE_ROOT
  // setup hists
  m_shape->setupHists();
#endif
  // set some constants
  m_rRange   = m_shape->extension();
  m_deltaR   = m_rRange / static_cast<double>(m_widthDiv);
  m_firstR   = m_deltaR/2. - m_rRange;

  m_phiRange = fabs(std::atan(m_rRange/rPerp));
  m_deltaPhi = m_phiRange / static_cast<double>(m_phiDiv);
  m_firstPhi = m_phiRange - m_deltaPhi/2.;

 #ifdef DEBUG
  std::cout << "[RadialProfile] constructed with:" << std::endl;
  std::cout << "                radius of detector:            " << rPerp << std::endl;
  std::cout << "                max radial cylinder extension: " << m_rRange << std::endl;
  std::cout << "                number of radial divisions:    " << m_widthDiv << std::endl;
  std::cout << "                radial step size:              " << m_deltaR << std::endl;
  std::cout << "                radial stepping starts at:     " << m_firstR << std::endl;
  std::cout << "                number of azimuthal divisions: " << m_phiDiv << std::endl;
  std::cout << "                azimuthal range:               " << m_phiRange*1000. << " x 10^-3" << std::endl;
  std::cout << "                azimuthal stepping size:       " << m_deltaPhi*1000. << " x 10^-3" << std::endl;
  std::cout << "                azimuthal stepping starts at:  " << m_firstPhi*1000. << " x 10^-3" << std::endl;
  
  m_shape->print();
 #endif
}
  
RadialProfile::~RadialProfile()
{
  if ( m_shape != 0 ) delete m_shape;
}

const container_t& RadialProfile::distribute(const particle_t& particle,
					     const SpacePoint& impact,
					     const SpacePoint& origin)
{
  return this->distribute(particle.e(),particle.pseudorapidity(),
			  particle.phi(),impact,origin);
}

const container_t& RadialProfile::distribute(double e,double eta,double phi,
					     const SpacePoint& impact,
					     const SpacePoint& /*origin*/)
{
  this->reset();
  // get distance of impact from vertex
  double rVtx(impact.distance());
  double rVtx2(rVtx*rVtx);
  // get angles
  double theta(std::acos(std::tanh(eta)));
  double firstPhi(phi-m_firstPhi);
  double lastPhi(phi+m_phiRange);

  // loops
  double r(m_firstR);
  double esum(0.);

  for ( ; r < m_rRange; r+=m_deltaR )
    {
      double da(std::atan(fabs(r)/rVtx));
      double atheta = r < 0. ? theta-da : theta+da;
      double sinTheta(std::sin(atheta));
      double cosTheta(std::cos(atheta));
      double aphi(firstPhi);
      double l(std::sqrt(rVtx2+r*r));
      //   std::cout << "[RadialProfile::distribute] (r_cyl,r_vtx,l,theta,phi) = (" 
      //		<< r << "," << rVtx << "," << l << "," << atheta << "," 
      //		<< aphi << ")" << std::endl;
      particle_t p;
      for ( ; aphi < lastPhi; aphi += m_deltaPhi )
	{
	  double x(l*sinTheta*std::cos(aphi));
	  double y(l*sinTheta*std::sin(aphi));
	  double z(l*cosTheta);
	  //SpacePoint t(x,y,z);
	  // check containment
	  double ep((m_shape->operator())(std::abs(r),e)*e); // access profile
    p.reset_momentum(ep*x/l,ep*y/l,ep*z/l,ep);
    m_profile.push_back(p);
	  //m_profile.push_back(particle_t(ep*x/l,ep*y/l,ep*z/l,ep));
	  esum += ep;
	}
    }
  // normalize integral energy
  double fe(e/esum);
  container_t::iterator fPart(m_profile.begin());
  container_t::iterator lPart(m_profile.end());
  //  esum = 0.;
  for ( ; fPart != lPart; ++fPart ) { (*fPart) *= fe; /*esum += (*fPart).e();*/ }
  return m_profile;
}

RadialProfile* RadialProfile::clone() const
{
  return new RadialProfile(*this);
}
#ifdef HAVE_ROOT
bool RadialProfile::writeHists()
{
  m_shape->writeHists();
}
#endif
//
// --------------------------------------------- Smearing
//
RadialSmearing::RadialSmearing()
  : m_emGrid(new Grid(Defaults::m_nEtaBinsEM,
		      Defaults::m_etaMinEM,
		      Defaults::m_etaMaxEM,
		      Defaults::m_nPhiBinsEM))
  , m_hadGrid(new Grid(Defaults::m_nEtaBinsHAD,
		       Defaults::m_etaMinHAD,
		       Defaults::m_etaMaxHAD,
		       Defaults::m_nPhiBinsHAD))
  , m_emProfile(0)
  , m_hadProfile(0)
  , m_detector(new CylindricalCalo(Defaults::m_caloRadius,
				   Defaults::m_caloHalfWidth,
				   Defaults::m_caloEtaMin,
				   Defaults::m_caloEtaMax))
  , m_oneGrid(false)
  , m_pseudoJets(0)
  , m_emPseudoJets(0)
  , m_hadPseudoJets(0)
{
  // get perpendicular radius
  CylindricalCalo* pCalo = dynamic_cast<CylindricalCalo*>(m_detector);
  double radius(0.);
  if ( pCalo != 0 )
    { 
      radius = pCalo->radius();  
 #ifdef DEBUG
      std::cout << "[RadialSmearing] use CylindricalCalo with radius " << radius << " mm" << std::endl;
 #endif
    }
  else
    {
      SpacePoint 
	p(m_detector->impact(particle_t(0.,10.,10.,std::sqrt(2.)*10.)));
      radius = std::sqrt(p.x()*p.x()+p.y()*p.y());
 #ifdef DEBUG
      std::cout << "[RadialSmearing] use detector inner radius " << radius << " mm" << std::endl;
 #endif
    }
  // new defaults
  //  m_emProfile = new RadialProfile(radius,30,32,RadialGaussShape(20.));
  //  m_hadProfile = new RadialProfile(radius,30,32,RadialGaussShape(40.));
#ifdef HAVE_ROOT
  m_emProfile = new RadialProfile(radius,30,32,RadialExpShapeEM(80.*Units::mm,"EMCalo"));
  m_hadProfile = new RadialProfile(radius,30,32,RadialExpShapeHAD(150.*Units::mm,"HADCalo"));
#else
  m_emProfile = new RadialProfile(radius,30,32,RadialExpShapeEM(80.*Units::mm));
  m_hadProfile = new RadialProfile(radius,30,32,RadialExpShapeHAD(150.*Units::mm));
#endif			  

  m_pseudoJets.reserve(100000);
  m_emPseudoJets.reserve(100000);
  m_hadPseudoJets.reserve(100000);

 #ifdef DEBUG
  std::cout << "[RadialSmearing] constructed with defaults! " << std::endl;
 #endif
}

RadialSmearing::RadialSmearing(const Grid& grid) 
  : m_emGrid(new Grid(grid))
  , m_hadGrid(new Grid(grid))
  , m_emProfile(0)
  , m_hadProfile(0)
  , m_detector(new CylindricalCalo(Defaults::m_caloRadius,
				   Defaults::m_caloHalfWidth,
				   Defaults::m_caloEtaMin,
				   Defaults::m_caloEtaMax))
  , m_oneGrid(false)
  , m_pseudoJets(0)
  , m_emPseudoJets(0)
  , m_hadPseudoJets(0)  
  
{ }

RadialSmearing::RadialSmearing(const Grid& emGrid,const Grid& hadGrid)
  : m_emGrid(new Grid(emGrid))
  , m_hadGrid(new Grid(hadGrid))
  , m_emProfile(0)
  , m_hadProfile(0)
  , m_detector(new CylindricalCalo(Defaults::m_caloRadius,
				   Defaults::m_caloHalfWidth,
				   Defaults::m_caloEtaMin,
				   Defaults::m_caloEtaMax))
  , m_oneGrid(false)
  , m_pseudoJets(0)
  , m_emPseudoJets(0)
  , m_hadPseudoJets(0)  
  
{ }


RadialSmearing::RadialSmearing(const RadialSmearing& module)
  : m_emGrid(module.m_emGrid->clone())
  , m_hadGrid(module.m_hadGrid->clone())
  , m_emProfile(module.m_emProfile->clone())
  , m_hadProfile(module.m_hadProfile->clone())
  , m_detector(module.m_detector->clone())
  , m_oneGrid(module.m_oneGrid)
  , m_pseudoJets(module.m_pseudoJets)
  , m_emPseudoJets(module.m_emPseudoJets)
  , m_hadPseudoJets(m_hadPseudoJets)
{ }

RadialSmearing::~RadialSmearing()
{
  if ( m_emGrid != 0 )     delete m_emGrid;
  if ( m_hadGrid != 0 )    delete m_hadGrid;
  if ( m_emProfile != 0 )  delete m_emProfile;
  if ( m_hadProfile != 0 ) delete m_hadProfile;
  if ( m_detector != 0 )   delete m_detector;
}

const container_t& RadialSmearing::smear(const particle_t& input,bool /*reset*/) {
  this->reset();
  
  if ( m_detector->inAcceptance(input) )
    {
      if ( this->isEM(input) ) 
        {
          this->fillEM(input);//,m_emGrid,m_emProfile,m_detector);
        }
      else
        {
          this->fillHAD(input);//,m_hadGrid,m_hadProfile,m_detector);
        }
    }
  
  // final copy
  m_emPseudoJets  = m_emGrid->allPseudoJets();
  m_hadPseudoJets = m_hadGrid->allPseudoJets(); 

  m_pseudoJets.insert(m_pseudoJets.end(),
		      this->emPseudoJets().begin(),
		      this->emPseudoJets().end());
  m_pseudoJets.insert(m_pseudoJets.end(),
		      this->hadPseudoJets().begin(),
		      this->hadPseudoJets().end());
  // return reference
  return m_pseudoJets;

}

const container_t& RadialSmearing::smear(const container_t& input,bool /*reset*/)
{
  // loop input particle
  this->reset();
  container_t::const_iterator fPart(input.begin());
  container_t::const_iterator lPart(input.end());
  for ( ; fPart != lPart; ++fPart )
    {
      if ( m_detector->inAcceptance(*fPart) )
	{
	  if ( this->isEM(*fPart) ) 
	    {
	      this->fillEM(*fPart);//,m_emGrid,m_emProfile,m_detector);
	    }
	  else
	    {
	      this->fillHAD(*fPart);//,m_hadGrid,m_hadProfile,m_detector);
	    }
	}
    }
  // final copy
  m_emPseudoJets  = m_emGrid->allPseudoJets();
  m_hadPseudoJets = m_hadGrid->allPseudoJets(); 

  m_pseudoJets.insert(m_pseudoJets.end(),
		      this->emPseudoJets().begin(),
		      this->emPseudoJets().end());
  m_pseudoJets.insert(m_pseudoJets.end(),
		      this->hadPseudoJets().begin(),
		      this->hadPseudoJets().end());
  // return reference
  return m_pseudoJets;
}

bool RadialSmearing::fillEM(const particle_t& particle)
{

  // get split energy
  //  SpacePoint p(m_detector->impact(particle));
  const container_t& sharing = 
    m_emProfile->distribute(particle,m_detector->impact(particle));

  // project split energies onto tower grid
  container_t::const_iterator fShare(sharing.begin());
  container_t::const_iterator lShare(sharing.end());
  for ( ; fShare != lShare; ++fShare ) 
    { 
      // energy smeared outside em acceptance is collected in HAD
      if ( !m_emGrid->add(*fShare) ) m_hadGrid->add(*fShare); 
    }

  // 
  return true;
}

bool RadialSmearing::fillHAD(const particle_t& particle)
{
  // get split energy
  const container_t& sharing = 
    m_hadProfile->distribute(particle,m_detector->impact(particle));

  // project split energies onto tower grid
  container_t::const_iterator fShare(sharing.begin());
  container_t::const_iterator lShare(sharing.end());
  for ( ; fShare != lShare; ++fShare ) { m_hadGrid->add(*fShare); }
  // 
  return true;
}

bool RadialSmearing::fill(const particle_t& particle,
			  Grid* grid,
			  RadialProfile* profile,
			  const DetectorDescription* detector)
{
  // check acceptance
  if ( !detector->inAcceptance(particle) ) return false;

  // get split energy
  SpacePoint p(detector->impact(particle));
  const container_t& sharing = profile->distribute(particle,p);

  // project split energies onto tower grid
  container_t::const_iterator fShare(sharing.begin());
  container_t::const_iterator lShare(sharing.end());
  for ( ; fShare != lShare; ++fShare ) { grid->add(*fShare); }

  // 
  return true;
}

bool RadialSmearing::finalize()
{
#ifdef HAVE_ROOT
  m_emProfile->writeHists();
  m_hadProfile->writeHists();
#endif
  return true;
}

double Math::poly(const std::vector<double>& pi,double x)
{
  // invalid input
  if ( pi.empty() ) return 0.;

  // calculations
  double eval(0.);
  double xp(1.);
  std::vector<double>::const_iterator fp(pi.begin());
  std::vector<double>::const_iterator lp(pi.end());
  for ( ; fp != lp; ++fp ) { eval += (*fp) * xp; xp *= x; }
  return eval;
}

double Math::poly(const std::vector<double>& pi,const std::vector<bool>& tags,
		  double x)
{
  // invalid input
  if ( pi.empty() ) return 0.;

  // calculations
  double eval(0.);
  double xp(1.);
  std::vector<double>::const_iterator fp(pi.begin());
  std::vector<double>::const_iterator lp(pi.end());
  size_t i(0);
  for ( ; fp != lp; ++fp, ++i ) 
    { 
      if ( tags[i] ) eval += (*fp) * xp; 
      xp *= x; 
    }
}

// #undef DEBUG
