// -*-c++-*-
#ifndef DETECTORMODEL_H_
#define DETECTORMODEL_H_

#include "fastjet/PseudoJet.hh"

#ifdef HAVE_ROOT
#include "TH2D.h"
#endif

#include <vector>
#include <cmath>
#include <limits>
#ifdef HAVE_ROOT
#include <string>
#include <map>
#endif

///////////////////////////////////////////////////////////////////////////////
//
// This code provides a crude lateral energy smearing for particle level MC.
//

namespace DetectorModel
{

  typedef fastjet::PseudoJet      particle_t;
  typedef std::vector<particle_t> container_t;
  typedef size_t                  index_t;

  struct Units
  {
    static double GeV;
    static double MeV;
    static double cm;
    static double mm;
  };

  struct Math
  {
    static double poly(const std::vector<double>& pi,double x);
    static double poly(const std::vector<double>& pi,
		       const std::vector<bool>& tags,
		       double x);
  };

  struct Defaults
  {
    static size_t m_nEtaBinsEM;
    static double m_etaMinEM;
    static double m_etaMaxEM;
    static size_t m_nPhiBinsEM;

    static size_t m_nEtaBinsHAD;
    static double m_etaMinHAD;
    static double m_etaMaxHAD;
    static size_t m_nPhiBinsHAD;

    static double m_caloRadius;
    static double m_caloHalfWidth;
    static double m_caloEtaMin;
    static double m_caloEtaMax;
  };

  class DetectorDescription;
  class Grid;
  class RadialProfile;

  ////////////////////////////////
  // Radial Energy Distribution //
  ////////////////////////////////

  class RadialSmearing
  {
  public:

    // construct with defaults
    RadialSmearing();
    // construct with client provided detector models
    RadialSmearing(const DetectorDescription& detector);
    // construct with client provided grid
    RadialSmearing(const Grid& grid);
    RadialSmearing(const Grid& emGrid,const Grid& hadGrid);
    // all freedom
    RadialSmearing(const DetectorDescription& detector,const Grid& emGrid,
		   const Grid& hadGrid);
    // copy constructor
    RadialSmearing(const RadialSmearing& module);

    virtual ~RadialSmearing();

    virtual const container_t& smear(const particle_t& input,bool reset=true);
    virtual const container_t& smear(const container_t& input,bool reset=true);

    virtual const container_t& allPseudoJets() const;
    virtual const container_t& emPseudoJets()  const;
    virtual const container_t& hadPseudoJets() const;

    virtual void setDetectorDescription(const DetectorDescription& detector);
    virtual void setEmGrid(const Grid& grid);
    virtual void setHadGrid(const Grid& grid);
    virtual void setEmProfile(const RadialProfile& profile);
    virtual void setHadProfile(const RadialProfile& profile);

    virtual RadialSmearing* clone() const;

    virtual bool finalize();

  protected:

    Grid* m_emGrid;
    Grid* m_hadGrid;
    RadialProfile* m_emProfile;
    RadialProfile* m_hadProfile;

    DetectorDescription* m_detector;
    bool m_oneGrid;

    // caching
    container_t m_pseudoJets;
    container_t m_emPseudoJets;
    container_t m_hadPseudoJets;

    void reset();
    bool fillEM(const particle_t& particle);
    bool fillHAD(const particle_t& particle);
    bool fill(const particle_t& particle,
	      Grid* grid,
	      RadialProfile* profile,
	      const DetectorDescription* detector);

    bool isEM(const particle_t& particle) const;

#ifdef HAVE_ROOT
    // debug plots
#endif
  };

  ////////////////
  // Tower Grid //
  ////////////////

  class Grid
  {
  public:

    // default grid
    Grid();
    Grid(const Grid& grid);
    // client provided grid
    Grid(size_t nEta,double etaMin,double etaMax,size_t nPhi);
    virtual ~Grid();

    virtual bool index(double eta,double phi,index_t& index) const;
    virtual bool index(index_t iEta,index_t iPhi,index_t& index) const;
    virtual bool phiIndex(double phi,index_t& index)         const;
    virtual bool phiIndex(index_t index,index_t& iPhi)   const;
    virtual bool etaIndex(double eta,index_t& index)         const;
    virtual bool etaIndex(index_t index,index_t& iEta)   const;

    virtual index_t nEta() const;
    virtual index_t nPhi() const;

    virtual double energy(index_t iEta,index_t iPhi) const;
    virtual double energy(index_t index)             const;
    virtual double energy(double eta,double phi)     const;

    virtual double eta(index_t index) const;
    virtual double phi(index_t index) const;

    virtual particle_t pseudoJet(index_t index)              const;
    virtual particle_t pseudoJet(index_t iEta,index_t iPhi)  const;
    virtual particle_t pseudoJet(double eta,double phi)      const;
    virtual particle_t pseudoJet(double e,double eta,double phi) const;

    virtual container_t allPseudoJets(bool zeroSuppression=true) const;

    virtual bool add(double e,double eta,double phi);
    virtual bool add(double px,double py,double pz,double e);
    virtual bool add(const particle_t& particle);

    virtual void reset();

    virtual Grid* clone() const;

  protected:

    std::vector<double> m_grid;

    size_t m_nBins;
    size_t m_nEta;
    size_t m_nPhi;
    double m_etaMin;
    double m_etaMax;
    double m_dEta;
    double m_dPhi;
    double m_twoPi;
    
    container_t fillNoZero() const;
    container_t fillAll()    const;
  };

  //////////////////////
  // Geometry Helpers //
  //////////////////////

  class SpacePoint;
  class Line;

  struct Tools
  {
    static double distance(const SpacePoint& pointA,const SpacePoint& pointB);
    static double distance(const Line& line,const SpacePoint& point);
    static double parallel(const Line& line,const SpacePoint& point);
    static double perp(const Line& line,const SpacePoint& point);
  };

  class SpacePoint
  {
  public:
    SpacePoint();
    SpacePoint(const SpacePoint& point);
    SpacePoint(double x,double y,double z);
    virtual ~SpacePoint();

    virtual double x() const;
    virtual double y() const;
    virtual double z() const;
    virtual double distance(const SpacePoint& point) const;
    virtual double distance() const;

  protected:
    double m_x;
    double m_y;
    double m_z;
    double m_rad;
  };

  class Line
  {
  public:
    Line();
    Line(const Line& line);
    Line(const SpacePoint& a,const SpacePoint& b);
    Line(const SpacePoint& a,double dirX,double dirY,double dirZ);
    virtual ~Line();

    virtual const SpacePoint& origin() const;
    virtual double length() const;

    virtual double distance(const SpacePoint& point) const;
    virtual double along(const SpacePoint& point) const;
    virtual double parallel(const SpacePoint& point) const;
    virtual double perp(const SpacePoint& point) const;

    virtual double dirX() const;
    virtual double dirY() const;
    virtual double dirZ() const;

  protected:
    SpacePoint m_origin;
    double     m_dirX;
    double     m_dirY;
    double     m_dirZ;
    double     m_length;
  };

  //////////////////////////
  // Detector Description //
  //////////////////////////

  class DetectorDescription
  {
  public:
    // default Geometry
    DetectorDescription();
    virtual ~DetectorDescription();

    virtual bool inAcceptance(double e,double eta,double phi,double m=0.) 
      const = 0;
    virtual bool inAcceptance(const particle_t& particle) const;

    virtual SpacePoint impact(double e,double eta,double phi,double m=0.) 
      const = 0;

    virtual DetectorDescription* clone() const = 0;


    virtual SpacePoint impact(const particle_t& particle) const;

  };

  class CylindricalCalo : public DetectorDescription
  { 
  public:
    CylindricalCalo(double r,double z,double etaMin,double etaMax);
    virtual ~CylindricalCalo();

    virtual bool inAcceptance(double e,double eta,double phi,double m) const;
    virtual SpacePoint impact(double e,double eta,double phi,double m) const;
    virtual DetectorDescription* clone() const;
 
    virtual double radius() const;
    virtual double dz() const;

  protected:
    double m_radius;
    double m_width;
    double m_etaMin;
    double m_etaMax;

    double m_etaBoundFwd;
    double m_etaBoundBwd;    

  };

  /////////////////////

  class RadialShape
  {
  public:
#ifndef HAVE_ROOT
    RadialShape(double width=160.*Units::mm);
#else
    RadialShape(double width=160*Units::mm,const std::string& name="NONE");
#endif
    virtual ~RadialShape();

    virtual double operator()(double r,double e) const = 0;
    virtual double evaluate(double r,double e=0.);

    virtual RadialShape* clone() const = 0;

    virtual double extension() const;

    virtual void print() const = 0;
#ifdef HAVE_ROOT
  public:
    virtual bool setupHists() = 0;
    virtual bool writeHists() = 0;
    virtual const std::string& name() const;
#endif
  protected:

    void setExtension(double width);
    double m_extension;
#ifdef HAVE_ROOT
  protected:
    bool m_histBooked;

  private:
    std::string m_name;
#endif
  };

  class RadialGaussShape : public RadialShape
  {
  public:
#ifndef HAVE_ROOT
    RadialGaussShape(double width=40.*Units::mm);
    RadialGaussShape(double width,double rMax);
#else
    RadialGaussShape(double width=40.*Units::mm,const std::string& name="Gauss");
    RadialGaussShape(double width,double rMax,const std::string& name="Gauss");
#endif
    RadialGaussShape(const RadialGaussShape& shape);
    virtual ~RadialGaussShape();

    virtual double operator()(double r,double e=0.) const;
    virtual RadialShape* clone() const;
    virtual void print() const;
#ifdef HAVE_ROOT
  public:
    virtual bool setupHists();
    virtual bool writeHists();
#endif
  protected:

    void setParameters();
    void normalize();

    double m_width;
    double m_width2;
    double m_const;

    static double m_probLvl;
#ifdef HAVE_ROOT
    TH2D* d_profile;
#endif
  };
  

			 /////////////////////////////////////////
  //////////////////////// Radial Electromagnetic Shower Shape //////////////////////////
  //                                                                                   //
  // This class implements a radial shower shape for electrons and photons. The para-  //
  // meterization is based on data and simulations from testbeam experiments for the   //
  // H1 liquid argon calorimeter. It should be considered typical for a precision EM   //
  // calorimeter, with no specific reference to the signal features of any other       //
  // specific one. The shape is best represented by the sum of two exponentials, one   //
  // describing the radial energy spread close to the shower axis (direction of flight //
  // of the incoming particle) and the other describing the wider spread from the      //
  // later stage of the EM shower. Usually these profiles are only slowly varying with //
  // the incoming particle energy, so here the lateral shower profile from 30 GeV      //
  // electrons is used (no energy dependence).                                         //
  //                                                                                   //
  ///////////////////////////////////////////////////////////////////////////////////////

  class RadialExpShapeEM: public RadialShape
  {
  public:
    // default constructor provides showers with a maximum radial extension of 80 mm
#ifndef HAVE_ROOT
    RadialExpShapeEM(double width=40.*Units::mm);
#else
    RadialExpShapeEM(double width=40.*Units::mm,const std::string& name="EMProfile");
#endif
    RadialExpShapeEM(const RadialExpShapeEM& shape);
    virtual ~RadialExpShapeEM();

    virtual double operator()(double r,double e=0.) const;
    virtual RadialShape* clone() const;
    virtual void print() const;
#ifdef HAVE_ROOT
  public:
    virtual bool setupHists();
    virtual bool writeHists();
#endif
  protected:
    
    void   normalize();
    double eval(double r)                    const;
    double exp0(double r)                    const;
    double exp1(double r)                    const;
    double integral(double rMin,double rMax) const;
    double integral(double width)            const;
    double containment(double fraction)      const;

    double m_norm;

    static double m_p0;
    static double m_p1;
    static double m_p2;
    static double m_p3;
#ifdef HAVE_ROOT
    TH2D* d_profile;
#endif
  };

  class RadialExpShapeHAD: public RadialShape
  {
  public:
#ifndef HAVE_ROOT
    RadialExpShapeHAD(double width=40.);
#else
    RadialExpShapeHAD(double width=40.,const std::string& name="HadProfile");
#endif
    RadialExpShapeHAD(const RadialExpShapeHAD& shape);
    virtual ~RadialExpShapeHAD();

    virtual double operator()(double r,double e) const;
    virtual double evaluate(double r,double e) const;
    virtual RadialShape* clone() const;
    virtual void print() const;
#ifdef HAVE_ROOT
  public:
    virtual bool setupHists();
    virtual bool writeHists();
#endif
  protected:

    // 
    void setParameters();
    //
    double refShape(double r)       const;
    double sDiff(double r)          const;
    double eParm(double e,double r) const;
    double rangeLimit(double e)     const;
    //
    double integral(double r,double e) const;
    double containment(double f,double e) const;

    // variables for energy extrapolation
    static double m_e0;      // energy range of reference data
    static double m_e1;
    static double m_slope;   // slope parameter of log approximation
    double m_kappa;          // equation constant
    // parametrization of shape difference
    static double m_r0;      // 1st radial range 0 - r0
    static double m_r1;      // 2nd (r0-r1) and third range (r1->infty)
    static double m_p00;     // 1st polynomial
    static double m_p01;     
    static double m_p02;
    static double m_p03;
    static double m_p10;     // 2nd polynominal
    static double m_p11;
    static double m_p12;
    static double m_p13;
    static double m_p20;     // 3rd polynomial
    static double m_p21;
    static double m_p22;
    // parameterization of reference shape
    static double m_refp0;
    static double m_refp1;
    static double m_refp2;
    static double m_refp3;
    // limitations of validity
    static double m_shapeLimitLo;
    static double m_shapeLimitHi;
    // parameterization of range limit
    static double m_rangep0;
    static double m_rangep1;
    static double m_rangep2;
    static double m_rangep3;
    static double m_minRange;

    // other constants
    double m_logE0;
    double m_logE1;
    //    double m_kappa;
    double m_appak;

    // containers for parameters
    std::vector<double> m_diff1;
    std::vector<double> m_diff2;
    std::vector<double> m_diff3;

    // range parametrization
    std::vector<double> m_range;
#ifdef HAVE_ROOT
    TH2D* d_profile;
    TH2D* d_range;
#endif
  };

  class RadialProfile
  {
  public:
    
    RadialProfile(double perp=1200.,
		  size_t widthDiv=30,
		  size_t phiDiv=64,
		  const RadialShape& distFct=RadialGaussShape(40.));
    virtual ~RadialProfile();

    virtual const container_t& distribute(const particle_t& particle,
					  const SpacePoint& impact,
					  const SpacePoint& origin=
					  SpacePoint());
    virtual const container_t& distribute(double e,
					  double eta,
					  double phi,
					  const SpacePoint& impact,
					  const SpacePoint& origin=
					  SpacePoint());
    virtual const container_t& last() const;
    virtual const RadialShape& shape() const;

    virtual void reset();

    virtual RadialProfile* clone() const;
#ifdef HAVE_ROOT
    virtual bool writeHists();
#endif

//   protected:

//     virtual Weights(double theta);
//     virtual double weightR(double r,double phi) const;
//     virtual double weightTheta(double theta,double phi) const;

  protected:
    
    container_t m_profile;
    RadialShape* m_shape;

    size_t m_widthDiv;
    size_t m_phiDiv;
    
    double m_rRange;
    double m_deltaR;
    double m_firstR;

    double m_phiRange;
    double m_deltaPhi;
    double m_firstPhi;

    //    std::vector<double> m_weights;
    //    std::vector<double> m_thetaBins;

    //    size_t index(double width,double phi) const;

    static double m_relPrec;

  };
}

//////////////////////
// Geometry Helpers //
//////////////////////

inline double 
DetectorModel::Tools::parallel(const Line& line,const SpacePoint& point)
{
  return 
    line.dirX()*(point.x()-line.origin().x()) +
    line.dirY()*(point.y()-line.origin().y()) +
    line.dirZ()*(point.z()-line.origin().z());
}

inline double
DetectorModel::Tools::perp(const Line& line,const SpacePoint& point)
{ return distance(line,point); }

////////////////
// SpacePoint //
////////////////

inline double DetectorModel::SpacePoint::x() const { return m_x; }
inline double DetectorModel::SpacePoint::y() const { return m_y; }
inline double DetectorModel::SpacePoint::z() const { return m_z; }
inline double 
DetectorModel::SpacePoint::distance(const SpacePoint& point) const
{ return Tools::distance(*this,point); }
inline double DetectorModel::SpacePoint::distance() const
{ return m_rad; }

//////////
// Line //
//////////

inline const DetectorModel::SpacePoint& DetectorModel::Line::origin() const
{ return m_origin; }

inline double DetectorModel::Line::length() const
{ return m_length; }

inline double DetectorModel::Line::dirX() const { return m_dirX; }
inline double DetectorModel::Line::dirY() const { return m_dirY; }
inline double DetectorModel::Line::dirZ() const { return m_dirZ; }

inline double DetectorModel::Line::along(const SpacePoint& point) const
{ return Tools::parallel(*this,point); }
inline double DetectorModel::Line::distance(const SpacePoint& point) const
{ return Tools::distance(*this,point); }

inline double DetectorModel::Line::perp(const SpacePoint& point) const
{ return this->distance(point); }
inline double DetectorModel::Line::parallel(const SpacePoint& point) const
{ return this->along(point); }

////////////////
// Tower Grid //
////////////////

inline bool 
DetectorModel::Grid::index(double eta,double phi,index_t& index) const
{
  index_t ii(0); index_t jj(0);
  return this->etaIndex(eta,ii) && this->phiIndex(phi,jj)
    ? this->index(ii,jj,index) : false;
}
inline bool 
DetectorModel::Grid::index(index_t iEta,index_t iPhi,index_t& index) const
{
  if ( iEta < m_nEta && iPhi < m_nPhi )
    {
      index = iEta * this->nPhi() + iPhi;
      return true;
    }
  return false;
}

inline DetectorModel::index_t DetectorModel::Grid::nEta() const 
{ return m_nEta; }
inline DetectorModel::index_t DetectorModel::Grid::nPhi() const 
{ return m_nPhi; }

inline bool DetectorModel::Grid::phiIndex(index_t index,index_t& iPhi) const
{ iPhi = index % m_nPhi; return iPhi < m_nPhi; }
inline bool DetectorModel::Grid::etaIndex(index_t index,index_t& iEta) const
{ iEta = index / m_nPhi; return iEta < m_nEta; }

inline double DetectorModel::Grid::energy(index_t iEta,index_t iPhi) const
{ 
  index_t index(0); 
  return this->index(iEta,iPhi,index) ? m_grid.at(index) : 0.; 
}
inline double DetectorModel::Grid::energy(index_t index)             const
{ return index < m_grid.size() ? m_grid.at(index) : 0.; }
inline double DetectorModel::Grid::energy(double eta,double phi)     const
{ 
  index_t index(0);
  return this->index(eta,phi,index) ? m_grid.at(index) : 0.;
}

inline double DetectorModel::Grid::eta(index_t index) const
{ 
  index_t ii(0); 
  return this->etaIndex(index,ii) 
    ? m_etaMin + (static_cast<double>(ii)+0.5) * m_dEta
    : std::numeric_limits<double>::infinity();
}  

inline double DetectorModel::Grid::phi(index_t index) const
{
  index_t jj(0);
  return this->phiIndex(index,jj)
    ? (static_cast<double>(jj)+0.5)*m_dPhi
    : 0.;
}

inline DetectorModel::particle_t 
DetectorModel::Grid::pseudoJet(index_t iEta,index_t iPhi) const
{ 
  index_t index(0); 
  return this->index(iEta,iPhi,index) 
    ? this->pseudoJet(index) 
    : particle_t(0.,0.,0.,0.);
}

inline DetectorModel::particle_t
DetectorModel::Grid::pseudoJet(double eta,double phi) const
{
  index_t index(0);
  return this->index(eta,phi,index)
    ? this->pseudoJet(index)
    : particle_t(0.,0.,0.,0.);
}

inline DetectorModel::particle_t 
DetectorModel::Grid::pseudoJet(double e,double eta,double phi) const
{
  return particle_t(e/std::cosh(eta)*std::cos(phi),
		    e/std::cosh(eta)*std::sin(phi),
		    e*std::tanh(eta),e);
}

inline DetectorModel::container_t 
DetectorModel::Grid::allPseudoJets(bool zeroSuppression) const
{
  return zeroSuppression ? this->fillNoZero() : this->fillAll();
}

inline void DetectorModel::Grid::reset() 
{ m_grid.clear(); m_grid.resize(m_nBins,0.); }

inline bool DetectorModel::Grid::add(const particle_t& particle) 
{
  return this->add(particle.e(),particle.pseudorapidity(),particle.phi());
}

inline bool DetectorModel::Grid::add(double e,double eta,double phi)
{
  index_t index(0);
  bool inGrid(this->index(eta,phi,index));
  if  ( inGrid ) m_grid[index] += e;
  return inGrid;
}

inline DetectorModel::Grid*
DetectorModel::Grid::clone() const
{
  return new Grid(*this);
}

//////////////////////
// Cylindrical Calo //
//////////////////////

inline double DetectorModel::CylindricalCalo::radius() const
{ return m_radius; }

inline double DetectorModel::CylindricalCalo::dz() const
{ return m_width; }

inline bool DetectorModel::CylindricalCalo::inAcceptance(double /*e*/,
							 double eta,
							 double /*phi*/,
							 double /*m*/) const
{
  return eta > m_etaMin && eta < m_etaMax;
}

///////////////////
// Radial Shapes //
///////////////////

inline double DetectorModel::RadialShape::extension() const
{ return m_extension; }
inline void DetectorModel::RadialShape::setExtension(double width)
{ m_extension = width; }

inline const DetectorModel::container_t& 
DetectorModel::RadialProfile::last() const
{ return m_profile; }

inline const DetectorModel::RadialShape& 
DetectorModel::RadialProfile::shape() const
{ return *m_shape; }

inline void DetectorModel::RadialProfile::reset()
{ m_profile.clear(); }

///////////////////////////
// Special Shower Shapes //
///////////////////////////

inline double DetectorModel::RadialExpShapeEM::eval(double r) const
{ return this->exp0(r) + this->exp1(r); }

inline 
double DetectorModel::RadialExpShapeEM::integral(double rMin,double rMax) const
{
  return 
    (this->exp0(rMax) - this->exp0(rMin))/m_p1 + 
    (this->exp1(rMax) - this->exp1(rMin))/m_p3;
}

inline
double DetectorModel::RadialExpShapeEM::integral(double width) const
{ return this->integral(0.,width); }


/////////////////////
// Radial Smearing //
/////////////////////

inline const DetectorModel::container_t&
DetectorModel::RadialSmearing::allPseudoJets() const
{ return m_pseudoJets; }

inline const DetectorModel::container_t&
DetectorModel::RadialSmearing::emPseudoJets() const
{ return m_emPseudoJets; }

inline const DetectorModel::container_t&
DetectorModel::RadialSmearing::hadPseudoJets() const
{ return m_hadPseudoJets; }

inline void DetectorModel::RadialSmearing::reset()
{
  m_pseudoJets.clear();
  m_emPseudoJets.clear();
  m_hadPseudoJets.clear();
  m_emGrid->reset();
  m_hadGrid->reset();
}

inline bool 
DetectorModel::RadialSmearing::isEM(const particle_t& particle) const
{
  int id(std::abs(particle.user_index()));
  return id == 11 || id == 22;
}

inline void DetectorModel::RadialSmearing::setEmGrid(const Grid& grid)
{ if ( m_emGrid != 0 ) delete m_emGrid; m_emGrid = grid.clone(); }

inline void DetectorModel::RadialSmearing::setHadGrid(const Grid& grid)
{ if ( m_hadGrid != 0 ) delete m_hadGrid; m_hadGrid = grid.clone(); }

inline void 
DetectorModel::RadialSmearing::setDetectorDescription(const 
						      DetectorDescription& 
						      detector)
{ if ( m_detector != 0 ) delete m_detector; m_detector = detector.clone(); }

inline void 
DetectorModel::RadialSmearing::setEmProfile(const RadialProfile& profile)
{ 
  if ( m_emProfile != 0 ) delete m_emProfile; 
  m_emProfile = profile.clone(); 
}

inline void 
DetectorModel::RadialSmearing::setHadProfile(const RadialProfile& profile)
{ 
  if ( m_hadProfile != 0 ) delete m_hadProfile; 
  m_hadProfile = profile.clone(); 
}

inline DetectorModel::RadialSmearing* 
DetectorModel::RadialSmearing::clone() const
{ return new RadialSmearing(*this); }

//////////////////
// Radial Shape //
//////////////////

#ifdef HAVE_ROOT
inline const std::string& DetectorModel::RadialShape::name() const { return m_name; }
#endif

#endif

