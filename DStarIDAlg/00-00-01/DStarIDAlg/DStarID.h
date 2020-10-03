#ifndef Physics_Analysis_DStarID_H
#define Physics_Analysis_DStarID_H 

#include <fstream>
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"

#include "VertexFit/IVertexDbSvc.h"
#include "AIDA/AIDA.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IAxis.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"
#include "AIDA/IHistogramFactory.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "EvtRecEvent/EvtRecEtaToGG.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "DTagTool/DTagTool.h"

//#include "Hist2Alg/ReadBeamInfFromDb.h"

#include "MeasuredEcmsSvc/IMeasuredEcmsSvc.h" // <<==== Ebeam

//-----------------------------------------
// EventNavigator
#include "EventNavigator/EventNavigator.h"

// MC data
#include "McTruth/McParticle.h"
//#include "McTruth/MdcMcHit.h" // Needed? Alexy has

// MDC reconstructed data
#include "MdcRecEvent/RecMdcKalTrack.h"
//#include "MdcRecEvent/RecMdcHit.h" // Needed? Alexy has

// EMC reconstructed data
#include "EmcRecEventModel/RecEmcShower.h"
//-----------------------------------------
#include <vector>
// Specify the namespace
using AIDA::IHistogram1D;
using AIDA::IHistogram2D;
using AIDA::IHistogram3D;
using AIDA::ICloud1D;
//#include "TH1F.h"

//MC stuff
#include "HepPDT/ParticleDataTable.hh"
//#include "PartPropSvc/PartPropSvc.h"
#include "GaudiKernel/IPartPropSvc.h"

class DStarID : public Algorithm {

public:
  DStarID(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

  void DumpTree();
  std::string ParticleName(int pdg);

private:
  int  m_irun;
  double  m_ebeam;
  int  m_myevt;
  IMeasuredEcmsSvc* ecmsSvc; // <<==== Ebeam
  
  HepPDT::ParticleDataTable* m_particleTable;
  // Functions
  //ReadBeamInfFromDb* m_readDb;
  HepLorentzVector myboost( const HepLorentzVector& p1,const HepLorentzVector& p2 );
  HepLorentzVector getP4(RecEmcShower* gTrk);
  HepLorentzVector gettrp4(RecMdcKalTrack* mdcKalTrack, int pid);
  HepLorentzVector gettrp4fac(RecMdcKalTrack* mdcKalTrack, int pid, double fac);
  int  MatchedSHID(int mcID,RecEmcShowerVector showers);
  bool  MatchedTK(int tkID, int mcPDG);
  int  MatchedTKID(int mcID,RecMdcKalTrackVector tracks);
  bool isTagTrue(DTagToolIterator iter_dtag);
  int  countPar(int parentID,int targetID);
  bool MCTKPID(int tkID, int mcPDG, int mcParPDG, int GParPDG );
  int MCTKPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG );
  bool MCPARPID(int tkID, int mcPDG, int mcParPDG             );
  bool MCSHPID(int tkID, int mcPDG, int mcParPDG, int GParPDG );
  int MCSHPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG );
  //
  bool Par2SH(int GParPDG, int mcPDG1, int mcPDG2,int tkID1, int tkID2);
  bool ETA3PI(int tkID1, int tkID2, int shID1, int shID2);
  bool ETAPGG(int tkID1, int tkID2, int shID1, int shID2);
  bool ETAP3P(int tkID1, int tkID2, int tkID3, int tkID4, int shID1, int shID2);
  bool RHOPP(int tkID1, int tkID2);
  bool ETAPRG(int tkID1, int tkID2, int shID1);
  bool KSPP(int tkID1, int tkID2);
  //-------------------------------------
  NTuple::Tuple* m_tuple1;
//
  NTuple::Item<int>     m_run;
  NTuple::Item<int>     m_EVT;

  NTuple::Item<int>     m_DSTp_piD0;
  NTuple::Item<int>     m_DSTm_piD0;
  NTuple::Item<int>     m_D0_kpi;
  NTuple::Item<int>     m_D0Bar_kpi;

  NTuple::Item<double>  m_DST_recMass;
  NTuple::Item<double>  m_DST_invMass;
  
  NTuple::Item<double>  m_D0_invMass;
  NTuple::Item<int>     m_D0_true;

  NTuple::Item<double>  m_slowPi_p;
  NTuple::Item<int>     m_slowPi_true;

  NTuple::Item<double>  m_HighestShw_Px;
  NTuple::Item<double>  m_HighestShw_Py;
  NTuple::Item<double>  m_HighestShw_Pz;

};


#endif 
