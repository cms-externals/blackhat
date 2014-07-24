#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "Interface/BH_Ampl.h"
#include "Interface/BH_interface.h"
#include <vector>


static BH::BH_interface *bhf=0;
static int bhlme_counter=0;

using namespace BH;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace OLP {
  class BH_Sherpa_Interface : public Virtual_ME2_Base {
    BH::BH_Ampl* p_BH;
  public:
    BH_Sherpa_Interface(const Process_Info& pi,const Flavour_Vector& flavs);
    ~BH_Sherpa_Interface();

    void Calc(const ATOOLS::Vec4D_Vector& mom);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);

    inline BH::BH_Ampl *Ampl() const { return p_BH; }
  };
}
using namespace OLP;

BH_Sherpa_Interface::BH_Sherpa_Interface(const Process_Info& pi, const Flavour_Vector& flavs) :
  Virtual_ME2_Base(pi, flavs)
{
  bhlme_counter++;

  std::vector<int> kflist;
  for (size_t i=0; i<flavs.size(); ++i) kflist.push_back(flavs[i].HepEvt());

  if (bhf==0) {
#ifndef BH_INTERFACE_BHSETTINGS
    bhf=new BH::BH_interface();
#else
    string bhfile("");
    Data_Reader reader(" ",";","!","=");
    if (reader.ReadFromFile(bhfile,"BH_SETTINGS_FILE")) {
      bhf=new BH::BH_interface(bhfile);
    }
    else bhf=new BH::BH_interface();
#endif
    rpa->gen.AddCitation(1,"The BlackHat library is described in \\cite{Berger:2008sj}.");
    bhf->set("Z_mass",Flavour(kf_Z).Mass());
    bhf->set("Z_width",Flavour(kf_Z).Width());
    bhf->set("W_mass",Flavour(kf_Wplus).Mass());
    bhf->set("W_width",Flavour(kf_Wplus).Width());
    double sin_th_2=MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"));
    if (MODEL::s_model->ScalarNumber(std::string("WidthScheme"))==1)
      sin_th_2=std::abs(MODEL::s_model->ComplexConstant(std::string("csin2_thetaW")));
    bhf->set("sin_th_2",sin_th_2);
    bhf->set("alpha_S",MODEL::s_model->ScalarFunction(std::string("alpha_S")));
    bhf->set("alpha_QED",MODEL::s_model->ScalarFunction(std::string("alpha_QED")));
  }
  
  if (pi.m_fi.m_sv=="FullColor")
    bhf->set("COLOR_MODE",std::string("full_color"));
  else if (pi.m_fi.m_sv=="LeadingColor")
    bhf->set("COLOR_MODE",std::string("leading_color"));
  else if (pi.m_fi.m_sv=="FullMinusLeadingColor")
    bhf->set("COLOR_MODE",std::string("full_minus_leading_color"));
  else if (pi.m_fi.m_sv!="")
    THROW(critical_error,"Invalid option '"+pi.m_fi.m_sv+"'");
  p_BH = bhf->new_ampl(kflist);

}

BH_Sherpa_Interface::~BH_Sherpa_Interface()
{
  delete p_BH;
  bhlme_counter--;
  if (bhlme_counter==0) {
    if (bhf) delete bhf;
    bhf=0;
  }
}

void BH_Sherpa_Interface::Calc(const Vec4D_Vector& mom)
{
  std::vector<std::vector<double> > moms(mom.size(), std::vector<double>(4, 0.0));
  for (size_t i=0; i<mom.size(); ++i) {
    for (size_t j=0; j<4; ++j) {
      moms[i][j]=mom[i][j];
    }
  }
  BH::BHinput input(moms, sqrt(m_mur2));
  bhf->operator()(input);

  m_res.Finite() = p_BH->get_finite();
  m_res.IR()     = p_BH->get_single_pole();
  m_res.IR2()    = p_BH->get_double_pole();
}

double BH_Sherpa_Interface::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{   
  return 4.*M_PI;
}

DECLARE_VIRTUALME2_GETTER(BH_Sherpa_Interface_Getter,"BH_Sherpa_Interface")
Virtual_ME2_Base *BH_Sherpa_Interface_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="BlackHat") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl=pi.ExtractFlavours();
    std::vector<int> kfvector;
    for (size_t i=0; i<fl.size(); ++i) kfvector.push_back(fl[i].HepEvt());
    BH_Sherpa_Interface *bh(new BH_Sherpa_Interface(pi, fl));
    if (bh->Ampl()->get_order_qed()!=pi.m_oew) return NULL;
    if (bh->Ampl()->get_order_qcd()+1!=pi.m_oqcd) return NULL;
    return bh;
  }
  return NULL;
}
