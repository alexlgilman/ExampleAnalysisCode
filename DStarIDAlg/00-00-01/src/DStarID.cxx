/* DStarIDAlg | 30/9/2020 | Author: Alex Gilman (Many functions written by Dr. Hajime Muramatsu)
 * This code identifies D^{*+} candidates through the decay chain D^{*+}-> \pi^{+} D^0, D^0->K^-\pi^+.
 * A single candidate is chosen based on the the invariant mass recoiling against the D^{*+} (recoil mass or mRec).
 * PDG codes are from https://pdg.lbl.gov/2006/reviews/pdf-files/montecarlo-web.pdf
 *
 * What's saved in the ntuple for each event:
 * > The measured E_CM, Beta Vector, and Run/Event No. 
 * > Truth-level information for simulated events
 * > Best D* Candidate  recoil mass and invariant mass
 * > Associated D0 Candidate 3-momentum and truth-matching
 * > Slow D* transition pion magnitude of 3-momentum, and truth-matching information of the associated 
 * > 3-momentum of the highest energy shower in the event
 */

//Dependencies

#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/AlgFactory.h"

#include "GaudiKernel/ISvcLocator.h"

#include "GaudiKernel/SmartDataPtr.h"

#include "GaudiKernel/IDataProviderSvc.h"

#include "GaudiKernel/IDataManagerSvc.h"

#include "GaudiKernel/PropertyMgr.h"

#include "SimplePIDSvc/ISimplePIDSvc.h"

#include "EventModel/EventModel.h"

#include "EventModel/EventHeader.h"

#include "EventModel/Event.h"

#include "TrigEvent/TrigEvent.h"

#include "TrigEvent/TrigData.h"

#include "EvTimeEvent/RecEsTime.h"

#include "EvtRecEvent/EvtRecEvent.h"

#include "EvtRecEvent/EvtRecTrack.h"

#include "DstEvent/TofHitStatus.h"

#include "EvtRecEvent/EvtRecVeeVertex.h"

#include "EvtRecEvent/EvtRecPi0.h"

#include "EvtRecEvent/EvtRecEtaToGG.h"

//-----------------------------------------
// MC stuff
#include "GaudiKernel/IPartPropSvc.h"

// EventNavigator
#include "EventNavigator/EventNavigator.h"

// MC data
#include "McTruth/McParticle.h"

// MDC reconstructed data
#include "MdcRecEvent/RecMdcKalTrack.h"

// EMC reconstructed data
#include "EmcRecEventModel/RecEmcShower.h"

// Muon
#include "MucRecEvent/RecMucTrack.h"
 //-----------------------------------------

#include "TMath.h"

#include "GaudiKernel/INTupleSvc.h"

#include "GaudiKernel/NTuple.h"

#include "GaudiKernel/Bootstrap.h"

#include "GaudiKernel/IHistogramSvc.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "CLHEP/Vector/TwoVector.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#include "CLHEP/Geometry/Point3D.h"

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D < double > HepPoint3D;
#endif

// Vertex (and K-FITS)
#include "VertexFit/KinematicFit.h"

#include "VertexFit/VertexFit.h"

#include "VertexFit/SecondVertexFit.h"

#include "VertexFit/Helix.h"

#include "VertexFit/IVertexDbSvc.h"

#include "VertexFit/WTrackParameter.h"

#include "ParticleID/ParticleID.h"

#include "MucRecEvent/RecMucTrack.h"

#include "VertexFit/KalmanKinematicFit.h"

// PID used in DTag
#include "SimplePIDSvc/ISimplePIDSvc.h"

#include "MeasuredEcmsSvc/IMeasuredEcmsSvc.h" // <<==== Run-dependent Ebeam for 4180 data
 //-----------------------------------------
#include "AIDA/AIDA.h"

#include "AIDA/IAnalysisFactory.h"

#include "AIDA/ITree.h"

#include "AIDA/IHistogram1D.h"

#include "AIDA/IHistogramFactory.h"

#include "AIDA/IAxis.h"

#include "AIDA/IHistogram2D.h"

#include "AIDA/IHistogram3D.h"

#include "AIDA/IHistogramFactory.h"

using AIDA::IHistogram1D;
using AIDA::IHistogram2D;
using AIDA::IHistogram3D;
using AIDA::ICloud1D;
//std::auto_ptr<AIDA::IAnalysisFactory>  af( AIDA_createAnalysisFactory() );
//std::auto_ptr<AIDA::ITreeFactory> tf( af -> createTreeFactory() ); 
//std::auto_ptr<AIDA::ITree> tree( tf -> create() ); 
//std::auto_ptr<AIDA::IHistogramFactory> hf( af->createHistogramFactory( *tree ) ); 
//-----------------------------------------
#include <iostream>

#include <fstream>
 //#include <iostream.h>
//#include <fstream.h>
//#include <iomanip.h>

using namespace std;
//-------------------------------------------

// needed to convert the Cell ID
#include "Identifier/Identifier.h"

#include "/afs/ihep.ac.cn/users/a/agilman/boss7.0.3p02/workarea/Physics/DStarIDAlg/00-00-01/DStarIDAlg/DStarID.h"

#include <vector>

// Define Constants
const double m_photonmass = 0.;
const double m_electmass = 0.000510998910;
const double m_muonmass = 0.105658367;
const double m_pionmass = 0.13957018;
const double m_pi0mass = 0.1349766;
const double m_kaonmass = 0.493677;
const double m_ksmass = 0.497614;
const double m_etamass = 0.547853;
const double m_omgmass = 0.78265;
const double m_protonmass = 0.938272013;
const double m_phimass = 1.019460;
const double m_D0mass = 1.86480;
const double m_dpmass = 1.86961;
const double m_dsmass = 1.96830;
const double m_dstzero = 2.00696;
const double m_dstplus = 2.01026;
const double m_dsstmass = 2.1121;
const double m_jpsimass = 3.096916;

// const double velc = 299.792458;   // tof path unit in mm
typedef std::vector < int > Vint;
typedef std::vector < HepLorentzVector > Vp4;

// counter
static long m_cout_all(0);

/////////////////////////////////////////////////////////////////////////////

DStarID::DStarID(const std::string & name, ISvcLocator * pSvcLocator):
  Algorithm(name, pSvcLocator) {
    //Declare the properties  

  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode DStarID::initialize() {
  cout << "Initialize" << endl;
  MsgStream logmess(msgSvc(), name());

  logmess << MSG::INFO << "in initialize()" << endmsg;

  StatusCode status;

  NTuplePtr nt1(ntupleSvc(), "FILE1/ups");
  if (nt1)
    m_tuple1 = nt1;
  else {
    m_tuple1 = ntupleSvc() -> book("FILE1/ups", CLID_ColumnWiseTuple, "My first example");

    if (m_tuple1) {
      status = m_tuple1 -> addItem("EVT", m_EVT);
      status = m_tuple1 -> addItem("run", m_run);

      status = m_tuple1 -> addItem("DSTp_piD0", m_DSTp_piD0);
      status = m_tuple1 -> addItem("DSTm_piD0", m_DSTm_piD0);
      status = m_tuple1 -> addItem("D0_kpi", m_D0_kpi);
      status = m_tuple1 -> addItem("D0Bar_kpi", m_D0Bar_kpi);

      status = m_tuple1 -> addItem("DST_recMass", m_DST_recMass);
      status = m_tuple1 -> addItem("DST_invMass", m_DST_invMass);

      status = m_tuple1 -> addItem("D0_invMass", m_D0_invMass);
      status = m_tuple1 -> addItem("D0_true", m_D0_true);

      status = m_tuple1 -> addItem("slowPi_p", m_slowPi_p);
      status = m_tuple1 -> addItem("slowPi_true", m_slowPi_true);

      status = m_tuple1 -> addItem("HighestShw_Px", m_HighestShw_Px);
      status = m_tuple1 -> addItem("HighestShw_Py", m_HighestShw_Py);
      status = m_tuple1 -> addItem("HighestShw_Pz", m_HighestShw_Pz);
    } else {
      logmess << MSG::ERROR << "    Cannot book N-nTuple:" << long(m_tuple1) << endmsg;

      return StatusCode::FAILURE;
    }

  }

  logmess << MSG::INFO << "successfully return from initialize()" << endmsg;

  return StatusCode::SUCCESS;
  // **********************************************************************

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode DStarID::execute() {
  m_myevt++;

  StatusCode sc = StatusCode::SUCCESS;
  //save the events passed selection to a new file
  setFilterPassed(false);

  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < Event::EventHeader > eventHeader(eventSvc(), "/Event/EventHeader");
  int runNo = eventHeader -> runNumber(); //Record run number (negative for MC events)
  int eventNo = eventHeader -> eventNumber(); // Record event number

  m_run = runNo;
  m_EVT = eventNo;

  m_ebeam = 4.17861 / 2; // Default Beam Energies
  if (fabs(runNo) >= 47543 && fabs(runNo) <= 48170) {
    m_ebeam = 4.18877 / 2.0;
  } //4190
  else if (fabs(runNo) >= 48172 && fabs(runNo) <= 48713) {
    m_ebeam = 4.19890 / 2.0;
  } //4200
  else if (fabs(runNo) >= 48714 && fabs(runNo) <= 49239) {
    m_ebeam = 4.20921 / 2.0;
  } //4210
  else if (fabs(runNo) >= 49270 && fabs(runNo) <= 49787) {
    m_ebeam = 4.21874 / 2.0;
  } //4220
  else if (fabs(runNo) >= 32239 && fabs(runNo) <= 32849) {
    m_ebeam = (4320.34 - 0.00287 * fabs(runNo)) * 0.5 / 1000.;
  } //4230
  else if (fabs(runNo) >= 32850 && fabs(runNo) <= 33484) {
    m_ebeam = 4.22554 / 2.0;
  } //4230

  Hep3Vector m_beta; //Default Beta Vector
  m_beta.setX(0.011);
  m_beta.setY(0);
  m_beta.setZ(0);

  Hep3Vector xorigin(0, 0, 0); // Read in the Interaction Point vertex from the database
  IVertexDbSvc * vtxsvc;
  Gaudi::svcLocator() -> service("VertexDbSvc", vtxsvc);
  if (vtxsvc -> isVertexValid()) {

    double * vertex = vtxsvc -> PrimaryVertex();
    xorigin.setX(vertex[0]);
    xorigin.setY(vertex[1]);
    xorigin.setZ(vertex[2]);
  }

  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < RecMucTrackCol > mucTracks(eventSvc(), "/Event/Recon/RecMucTrackCol");
  SmartDataPtr < EvtRecVeeVertexCol > evtRecVeeVertexCol(eventSvc(), "/Event/EvtRec/EvtRecVeeVertexCol");
  SmartDataPtr < EvtRecPi0Col > recPi0Col(eventSvc(), "/Event/EvtRec/EvtRecPi0Col");
  SmartDataPtr < EvtRecEtaToGGCol > recEtaCol(eventSvc(), "/Event/EvtRec/EvtRecEtaToGGCol");

  ////////////////////////// Generator Truth Information ////////////////////////////////////
  /*
   * Loop over the generated decay tree to find the processes of interest: DST->Pi D0 and D0->KPi.
   */
  if (runNo < 0) {
    SmartDataPtr < Event::McParticleCol > mcParticleCol(eventSvc(), EventModel::MC::McParticleCol);
    if (mcParticleCol) {
      Event::McParticleCol::iterator iter_mc = mcParticleCol -> begin();
      for (; iter_mc != mcParticleCol -> end(); iter_mc++) {

        if (( * iter_mc) -> particleProperty() == 413) //DSTp
        {
          const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
          if (gc.size() == 2) {
            for (unsigned int D0 = 0; D0 < gc.size(); D0++) {
              if (gc[D0] -> particleProperty() == 421) //D0
              {
                for (unsigned int pi = 0; pi < gc.size(); pi++) {
                  if (D0 == pi) {continue;}
                  if (gc[pi] -> particleProperty() == 211) //Pi+
                  {
                    m_DSTp_piD0 = 1;
                  } //pi found
                } //pi loop
              } //D0 found
            } //D0 loop
          } // if 2 daughters
        } //if DSTp

        if (( * iter_mc) -> particleProperty() == -413) //DSTm
        {
          const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
          if (gc.size() == 2) {
            for (unsigned int D0 = 0; D0 < gc.size(); D0++) {
              if (gc[D0] -> particleProperty() == -421) //D0Bar
              {
                for (unsigned int pi = 0; pi < gc.size(); pi++) {
                  if (D0 == pi) {continue;}
                  if (gc[pi] -> particleProperty() == -211) //Pi-
                  {
                    m_DSTm_piD0 = 1;
                  } //pi found
                } //pi loop
              } //D0 found
            } //D0 loop
          } // if 2 daughters
        } //if DSTm

        if (( * iter_mc) -> particleProperty() == 421) //D0
        {
          const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
          if (gc.size() == 2) {
            for (unsigned int k = 0; k < gc.size(); k++) {
              if (gc[k] -> particleProperty() == -321) //K-
              {
                for (unsigned int pi = 0; pi < gc.size(); pi++) {
                  if (k == pi) {
                    continue;
                  }
                  if (gc[pi] -> particleProperty() == 211) //Pi+
                  {
                    m_D0_kpi = 1;
                  } //pi found
                } //pi loop
              } //k found
            } //k loop
          } // if 2 daughters
        } //if D0

        if (( * iter_mc) -> particleProperty() == -421) //D0Bar
        {
          const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
          if (gc.size() == 2) {
            for (unsigned int k = 0; k < gc.size(); k++) {
              if (gc[k] -> particleProperty() == 321) //K+
              {
                for (unsigned int pi = 0; pi < gc.size(); pi++) {
                  if (k == pi) {
                    continue;
                  }
                  if (gc[pi] -> particleProperty() == -211) //Pi-
                  {
                    m_D0Bar_kpi = 1;
                  } //pi found
                } //pi loop
              } //k found
            } //k loop
          } // if 2 daughters
        } //if D0Bar

      } //generated particle loop
    } //mc particle col
  } //if MC
  /////////////////////// Fill Defaults ////////////////////////////

  HepLorentzVector m_Ecm = 2 * m_ebeam * HepLorentzVector(0, 0, 0, 1); // in the cm frame

  double min_delta_mRec = 10000;

  m_DST_recMass = 0;
  m_DST_invMass = 0;

  m_D0_invMass = 0;
  m_D0_true = 0;

  m_slowPi_p = 0;
  m_slowPi_true = 0;

  m_HighestShw_Px = -5;
  m_HighestShw_Py = -5;
  m_HighestShw_Pz = -5;

  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Begin DTag

  DTagTool dtagTool;
  if (dtagTool.isDTagListEmpty()) {}
  dtagTool.setPID(true);
  DTagToolIterator iter_begin = dtagTool.modes_begin();
  DTagToolIterator iter_end = dtagTool.modes_end();
  vector < int > d0itindex = dtagTool.D0modes();

  for (int o = 0; o < d0itindex.size(); o++) // Loop over D0 Candidates
  {
    DTagToolIterator iter_dtag = dtagTool.modes_begin() + d0itindex[o];
    int m_mode = ( * iter_dtag) -> decayMode();
    if (m_mode != 0) {
      continue;
    } // Require D0->KPi Candidate
    int m_charm = ( * iter_dtag) -> charm();
    // Default for charged kaons and pions are just be "good tracks".
    //     // If PID is needed, require pdg=1 (=EvtRecDTag::Tight)
    int m_type = ( * iter_dtag) -> type();
    if (m_type != 1) {
      continue;
    } //Require PID be applied to Kaon/Pion daughters of the D candidate
    //cout << "PID" << endl;  
    if (( * iter_dtag) -> type() != EvtRecDTag::Tight) {
      continue;
    }
    //        //cout << "EvtRecDTag" << endl;    

    HepLorentzVector D0_p4 = ( * iter_dtag) -> p4();

    SmartRefVector < EvtRecTrack > tracks = ( * iter_dtag) -> tracks(); //List of tracks used in the D0 reconstruction
    SmartRefVector < EvtRecTrack > showers = ( * iter_dtag) -> showers(); //List of showers used in the D0 reconstruction
    SmartRefVector < EvtRecTrack > othertracks = ( * iter_dtag) -> otherTracks(); //List of tracks unused in the D0 reconstruction
    SmartRefVector < EvtRecTrack > othershowers = ( * iter_dtag) -> otherShowers(); //List of showers unused in the D0 reconstruction

    // D0 Candidate Tracks
    RecMdcKalTrack * k_mdcTrk = (tracks[0]) -> mdcKalTrack();
    k_mdcTrk -> setPidType(RecMdcKalTrack::kaon);
    RecMdcKalTrack * pi_mdcTrk = (tracks[1]) -> mdcKalTrack();
    pi_mdcTrk -> setPidType(RecMdcKalTrack::pion);

    // Truth Matching: Compare simulated detector hits with generated trajectories to match a track with a generated particle
    // The MCTKPIDCHG function can return the PDG code associated with a track or the code associated its parent. It can also be used to require that the PDG codes of the particle (2nd argument) and its parent (3rd argument) match certain values. See the MCTKPIDCHG function for more details.
    int k_pdg = 0;
    int k_parent_pdg = 0;
    int pi_pdg = 0;
    int pi_parent_pdg = 0;
    if (runNo < 0) {
      k_pdg = MCTKPIDCHG(k_mdcTrk -> trackId(), -1, 0, 0);
      k_parent_pdg = MCTKPIDCHG(k_mdcTrk -> trackId(), 0, -1, 0);
      pi_pdg = MCTKPIDCHG(pi_mdcTrk -> trackId(), -1, 0, 0);
      pi_parent_pdg = MCTKPIDCHG(pi_mdcTrk -> trackId(), 0, -1, 0);
    }

    ///////////////Find Slow Transition Pion for Constrctuing D* /////////////////////

    HepLorentzVector this_slowPi_p4;
    int this_slowPi_pdg = 0;
    int this_slowPi_parent_pdg = 0;

    double this_mRec = 0;
    double this_delta_mRec = 10000;
    for (int i = 0; i < othertracks.size(); i++) //loop through all tracks unused in the reconstruction of the D0 candidate
    {
      RecMdcKalTrack * slowPi_mdcTrk = (othertracks[i]) -> mdcKalTrack();
      slowPi_mdcTrk -> setPidType(RecMdcKalTrack::pion);
      if (slowPi_mdcTrk -> charge() == 0) {continue;}

      // Find the track vertex
      HepVector slowPi_a = slowPi_mdcTrk -> helix();
      HepSymMatrix slowPi_Ea = slowPi_mdcTrk -> err();
      HepPoint3D point0(0., 0., 0.);
      VFHelix slowPi_helixip(point0, slowPi_a, slowPi_Ea);
      HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
      slowPi_helixip.pivot(IP);
      HepVector slowPi_vecipa = slowPi_helixip.a();
      double slowPi_vz = slowPi_vecipa[3];
      double slowPi_vr = fabs(slowPi_vecipa[0]);

      //Get the reconstructed 4-momentum (2nd argument relates to mass hypothesis
      HepLorentzVector slowPi_p4 = gettrp4(slowPi_mdcTrk, 2);//

      // Apply typical "good track" requirements: 
      // > Track vertex is within 1 (10) cm perpendicular (parallel) to the beam
      if (slowPi_vr >= 1 || fabs(slowPi_vz) >= 10) {continue;}
      // > Require the tracked angle w.r.t. the beampipe satisfies |cos theta|<0.93
      if (fabs(cos(slowPi_mdcTrk -> theta())) >= 0.93) {continue;}
      double slowPi_charge = slowPi_mdcTrk -> charge();

      // Apply standard pion identification using the ParticleID package. This uses dEdx information from the MDC and corrected measurements from the TOF. No information from the EMC is used.
      //
      ParticleID * slowPi_kpi_pid = ParticleID::instance();
      slowPi_kpi_pid -> init();
      slowPi_kpi_pid -> setMethod(slowPi_kpi_pid -> methodProbability());
      slowPi_kpi_pid -> setChiMinCut(4);
      slowPi_kpi_pid -> setRecTrack(othertracks[i]);
      slowPi_kpi_pid -> usePidSys(slowPi_kpi_pid -> useDedx() | slowPi_kpi_pid -> useTofCorr());
      slowPi_kpi_pid -> identify(slowPi_kpi_pid -> onlyPion() | slowPi_kpi_pid -> onlyKaon());
      slowPi_kpi_pid -> calculate();

      double slowPi_probP = slowPi_kpi_pid -> probPion();
      double slowPi_probK = slowPi_kpi_pid -> probKaon();

      if (!(slowPi_probP > slowPi_probK && slowPi_probP > 0)) {continue;} // Standard pion ID

      // Truth Match the Slow Pion Candidate

      // Calculate the momentum of the slow pion candidate in the D* rest frame
      HepLorentzVector DST_p4 = slowPi_p4 + D0_p4;
      HepLorentzVector slowPi_p4_BOOSTED = myboost(slowPi_p4, DST_p4);

      // Cut on the boosted pion momentum (<60 MeV/c) and require the pion have opposite charge as the D0 kaon daughter to satisfy the decay chain hypothesis
      if (slowPi_p4_BOOSTED.vect().mag() >= 0.060 || slowPi_charge == k_mdcTrk -> charge()) {continue;}

      // Boost the D* momentum to the COM frame to calculate the recoil mass
      DST_p4.boost(-m_beta);
      HepLorentzVector con_tag_p4 = HepLorentzVector(DST_p4.vect(), sqrt(DST_p4.vect().mag2() + m_dstplus * m_dstplus));
      double rec_mass = (m_Ecm - con_tag_p4).mag(); //recoil mass

      // Pick a best candidate based on recoil mass
      // DISCLAIMER: Picking a single candidate between a D*D* and D*D hypothesis is likely not an optimal analysis strategy, but I'm doing it here for simplicity.
      double delta_mRec_DSTDST = rec_mass - m_dstplus;
      double delta_mRec_DSTD = rec_mass - m_dpmass;

      //if this D* (pi D0) candidate is closer in recoil mass to either a D*D or D*D* hypothesis, replace it as the best candidate
      if (min(fabs(delta_mRec_DSTDST), fabs(delta_mRec_DSTD)) < min_delta_mRec) {
        // Populate the variables of interest for this best candidate
        this_delta_mRec = min(fabs(delta_mRec_DSTDST), fabs(delta_mRec_DSTD));
        this_mRec = rec_mass;
        this_slowPi_p4 = slowPi_p4;
        if (runNo < 0) // Only for MC
        {
          this_slowPi_pdg = MCTKPIDCHG((othertracks[i]) -> trackId(), -1, 0, 0);
          this_slowPi_parent_pdg = MCTKPIDCHG((othertracks[i]) -> trackId(), 0, -1, 0);
        }

        //Set the best candidate recoil mass to that of this candidate	
        min_delta_mRec = this_delta_mRec;
      }
    } //slow pion loop

    if (this_delta_mRec == 10000) {continue;} // No best candidate found for this D0 candidate, so skip the rest

    // Find the Highest Energy Shower in the Event unused in the tag reconstruction
    /* NOTE: Since this code only examines the D0->K-Pi+ mode, this will be the same for every candidate, and so it would be more efficient to run over all unmatched showers outside this loop.
     * However, I've left it here, because this would not be the case for tag candidates with neutral daughters.
     */
    double highest_ShwEn = 0;
    double this_ShwPx = 0;
    double this_ShwPy = 0;
    double this_ShwPz = 0;

    for (int i = 0; i < othershowers.size(); i++) {
      if (!(othershowers[i]) -> isEmcShowerValid()) {
        continue;
      }
      RecEmcShower * emcShw = (othershowers[i]) -> emcShower();
      // Apply "good" shower requirements:
      // > In the "Barrel" (|cos theta| <0.80) region of the EMC and has a minimum energy of 25 MeV
      // OR
      // > In the "Good Endcap" (0.93>|cos theta| >0.86) region of the EMC and has a minimum energy of 50 MeV
      if (!((fabs(cos(emcShw -> theta())) < 0.8 && emcShw -> energy() > 0.025) || (fabs(cos(emcShw -> theta())) > 0.86 && fabs(cos(emcShw -> theta())) < 0.93 && emcShw -> energy() > 0.050))) {continue; }

      //Get the shower 4-momentum
      HepLorentzVector shw_p4 = getP4(emcShw);

      //Check if it's the highest energy, and replace variables of interest if so.
      if (shw_p4.vect().mag() > highest_ShwEn) {
        this_ShwPx = shw_p4.x();
        this_ShwPy = shw_p4.y();
        this_ShwPz = shw_p4.z();

        highest_ShwEn=shw_p4.vect().mag();
      }

    } //for all unmatched showers unused in the D0 candidate reconstruction

    //Populate the ntuple entries, since we have already required that this be the best D* candidate  

    m_DST_recMass = this_mRec;
    m_DST_invMass = (this_slowPi_p4 + D0_p4).mag();

    m_D0_invMass = D0_p4.mag();
    m_D0_true = 0;
    if ((m_D0_kpi && k_pdg == -321 && k_parent_pdg == 421 && pi_pdg == 211 && pi_parent_pdg == 421) || (m_D0Bar_kpi && k_pdg == 321 && k_parent_pdg == -421 && pi_pdg == -211 && pi_parent_pdg == -421)) {
      m_D0_true = 1;
    }

    m_slowPi_p = this_slowPi_p4.vect().mag();
    m_slowPi_true = 0;
    if ((this_slowPi_pdg == 211 && this_slowPi_parent_pdg == 413) || (this_slowPi_pdg == -211 && this_slowPi_parent_pdg == -413)) {
      m_slowPi_true = 1;
    }

    m_HighestShw_Px = this_ShwPx;
    m_HighestShw_Py = this_ShwPy;
    m_HighestShw_Pz = this_ShwPz;
  } //D0 candidate loop
  /////////////////////////////////////////

  //Write to the mtuple: 1 entry per event
  m_tuple1 -> write();
  return sc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode DStarID::finalize() {

  MsgStream logmess(msgSvc(), name());
  cout << "_/_/_/_/_/_/_/_/_/ all event- " << m_cout_all << endl;
  return StatusCode::SUCCESS;
}
// ************************** FUNCTIONS *****************************************************
HepLorentzVector DStarID::myboost(const HepLorentzVector & p1,
  const HepLorentzVector & p2) // boost particle p1 to the rest frame of p2 (boosted p1 returned as the function value)
{
  double am2(p2.m());
  double e((p1.t() * p2.t() - (p1.x() * p2.x() + p1.y() * p2.y() +
    p1.z() * p2.z())) / am2);
  double r((p1.t() + e) / (p2.t() + am2));
  return HepLorentzVector(p1.x() - r * p2.x(),
    p1.y() - r * p2.y(),
    p1.z() - r * p2.z(),
    e);
}

HepLorentzVector DStarID::getP4(RecEmcShower * gTrk) {

  double eraw = gTrk -> energy();
  double phi = gTrk -> phi();
  double the = gTrk -> theta();

  return HepLorentzVector(eraw * sin(the) * cos(phi),
    eraw * sin(the) * sin(phi),
    eraw * cos(the),
    eraw);
}
HepLorentzVector DStarID::gettrp4(RecMdcKalTrack * mdcKalTrack, int pid) {
  HepVector zhelix;
  double mass = 0;

  if (pid == 0) {
    zhelix = mdcKalTrack -> getZHelixE();
    mass = 0.000511;
  } else if (pid == 1) {
    zhelix = mdcKalTrack -> getZHelixMu();
    mass = 0.105658;
  } else if (pid == 2) {
    zhelix = mdcKalTrack -> getZHelix();
    mass = 0.139570;
  } else if (pid == 3) {
    zhelix = mdcKalTrack -> getZHelixK();
    mass = 0.493677;
  } else {
    zhelix = mdcKalTrack -> getZHelixP();
    mass = 0.938272;
  }

  double dr(0), phi0(0), kappa(0), dz(0), tanl(0);
  dr = zhelix[0];
  phi0 = zhelix[1];
  kappa = zhelix[2];
  dz = zhelix[3];
  tanl = zhelix[4];

  int charge = 0;
  if (kappa > 0.0000000001)
    charge = 1;
  else if (kappa < -0.0000000001)
    charge = -1;

  double pxy = 0;
  if (kappa != 0) pxy = 1.0 / fabs(kappa);

  double px = pxy * (-sin(phi0));
  double py = pxy * cos(phi0);
  double pz = pxy * tanl;

  double e = sqrt(pxy * pxy + pz * pz + mass * mass);

  return HepLorentzVector(px, py, pz, e);
}
HepLorentzVector DStarID::gettrp4fac(RecMdcKalTrack * mdcKalTrack, int pid, double fac) {
  HepVector zhelix;
  double mass = 0;

  if (pid == 0) {
    zhelix = mdcKalTrack -> getZHelixE();
    mass = 0.000511;
  } else if (pid == 1) {
    zhelix = mdcKalTrack -> getZHelixMu();
    mass = 0.105658;
  } else if (pid == 2) {
    zhelix = mdcKalTrack -> getZHelix();
    mass = 0.139570;
  } else if (pid == 3) {
    zhelix = mdcKalTrack -> getZHelixK();
    mass = 0.493677;
  } else {
    zhelix = mdcKalTrack -> getZHelixP();
    mass = 0.938272;
  }

  double dr(0), phi0(0), kappa(0), dz(0), tanl(0);
  dr = zhelix[0];
  phi0 = zhelix[1];
  kappa = zhelix[2];
  dz = zhelix[3];
  tanl = zhelix[4];

  int charge = 0;
  if (kappa > 0.0000000001)
    charge = 1;
  else if (kappa < -0.0000000001)
    charge = -1;

  double pxy = 0;
  if (kappa != 0) pxy = 1.0 / fabs(kappa);

  double px = pxy * (-sin(phi0));
  double py = pxy * cos(phi0);
  double pz = pxy * tanl;

  double e = sqrt(pxy * pxy + pz * pz + mass * mass);

  px = fac * px;
  py = fac * py;
  pz = fac * pz;
  e = sqrt(px * px + py * py + pz * pz + mass * mass);

  return HepLorentzVector(px, py, pz, e);
}
//int DStarID::MCGENQ(int parID, int children){
//  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), EventModel::MC::McParticleCol);
//  int numChildren=0;
//  if(mcParticleCol)
//  {
//	Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
//  SmartRefVector<Event::McParticle> gc;
//  for ( ; iter_mc != mcParticleCol->end(); iter_mc++) 
//        {
//          if((*iter_mc)->particleProperty() == parID)
//          {
//          gc=(*iter_mc)->daughterList();
//               {
//                  for(unsigned int child = 0; child < gc.size(); child++) 
//                  {
//	                  if(!(gc[child]->decayInFlight())){numChildren++}
//                    if(gc[trk1]->particleProperty()==daughterList[0])
//                    {
//                     for(unsigned int k1 = 0; k1 < gc.size(); k1++)
//		     {
//	               if(k1==pi){continue;} 
//                     if(gc[k1]->particleProperty()==321)
//		     {
//                     for(unsigned int k2 = 0; k2 < gc.size(); k2++)
//		     {
//	               if(k2==k1){continue;} 
//	               if(k2==pi){continue;} 
//                       if(gc[k2]->particleProperty()==-321)
//		       {
//			       m_Dsp_kkp=1;
//		       }//if second kaon found
//		     }//second kaon loop
//		     }//if first kaon found
//		     }//first kaon loop
//		  }//if pion found
//	       }//pion loop
//	       }//if 3 trax
//
//}

int DStarID::MCTKPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  int ismatched = 0;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcParPDG == -1 && (GParPDG == 0 || GParPDG == thebestGParent)) {
        return thebestParent;
      }
      if (mcParPDG == -1 && GParPDG == 0) {
        return thebestParent;
      }
      if (mcPDG == -1) {
        return thebestmatchedPDG;
      }

      if (mcPDG == 0 && mcParPDG == 0) {
        if (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
          thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG) {
          return 1;
        }
      }
      if (GParPDG == 0 && mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG) {
          return 1;
        }
      }
      if (mcPDG == 0 && GParPDG == 0) {
        if (thebestParent == mcParPDG) {
          return 1;
        }
      }
      if (mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          return 1;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = 1;
        }
      }
    }
  }

  return ismatched;
}
bool DStarID::MCPARPID(int tkID, int mcPDG, int mcParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  bool ismatched = false;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcPDG == 0) {
        if (thebestParent == mcParPDG) {
          ismatched = true;
        }
      } else if (mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG) {
          ismatched = true;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG) {
          ismatched = true;
        }
      }
    }
  }

  return ismatched;
}
bool DStarID::MCTKPID(int tkID, int mcPDG, int mcParPDG, int GParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  bool ismatched = false;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      int thebestGGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        if (!particles[i] -> decayFromGenerator()) {
          continue;
        }
        if (particles[i] -> decayInFlight()) {
          continue;
        }
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcParPDG < 0) {
        if (thebestmatchedPDG == mcPDG) {
          ismatched = true;
        }
      } else if (mcParPDG == 0) {
        if (fabs(thebestmatchedPDG) == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG ||
            thebestGGGGGGGGParent == GParPDG)) {
          ismatched = true;
        }
      } else {
        if (fabs(thebestmatchedPDG) == mcPDG && fabs(thebestParent) == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG ||
            thebestGGGGGGGParent == GParPDG)) {
          ismatched = true;
        }
      }
    }
  }

  return ismatched;
}

bool DStarID::MCSHPID(int tkID, int mcPDG, int mcParPDG, int GParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  bool ismatched = false;
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == tkID) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcParPDG == 0) {
        if (fabs(thebestmatchedPDG) == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          //if (temp_hits>20) {ismatched = true;}
          ismatched = true;
        }
      } else {
        if (fabs(thebestmatchedPDG) == mcPDG && fabs(thebestParent) == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = true;
        }
      }
    }
  }

  return ismatched;
}
int DStarID::MCSHPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  int ismatched = -1;
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == tkID) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }

      if (mcParPDG == -1) {
        return thebestParent;
      }
      if (mcPDG == -1) {
        return thebestmatchedPDG;
      }

      if (mcPDG == 0 && mcParPDG == 0) {
        if (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
          thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG) {
          return 1;
        }
      }
      if (GParPDG == 0 && mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG) {
          return 1;
        }
      }
      if (mcPDG == 0 && GParPDG == 0) {
        if (thebestParent == mcParPDG) {
          return 1;
        }
      }
      if (mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          return 1;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = 1;
        }
      }
    }
  }
  return ismatched;
}
bool DStarID::Par2SH(int GParPDG, int mcPDG1, int mcPDG2, int tkID1, int tkID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      int matchedPDG = -2;
      int Parent = -2;
      int ParentMCID = -2;
      int MC_sh1 = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          matchedPDG = particles[i] -> particleProperty();
          Parent = particles[i] -> mother().particleProperty();
          ParentMCID = particles[i] -> mother().trackIndex();
          MC_sh1 = particles[i] -> trackIndex();
        }
      }
      if (matchedPDG == mcPDG1 && Parent == GParPDG) {
        for (iter_mc = mcParticles -> begin(); iter_mc != mcParticles -> end(); iter_mc++) {
          if (ParentMCID == ( * iter_mc) -> trackIndex()) {
            const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
            if (gc.size() > 0) {
              for (unsigned int ii = 0; ii < gc.size(); ii++) {
                if (gc[ii] -> particleProperty() == mcPDG2 &&
                  MC_sh1 != gc[ii] -> trackIndex()) {
                  int MCtkID2 = -2;
                  RecEmcShowerVector tracks2 = navigator -> getEmcRecShowers(gc[ii]);
                  int temp_hits2 = 0;
                  for (unsigned int z = 0; z < tracks2.size(); z++) {
                    int tkHITS = navigator -> getMcParticleRelevance(tracks2[z], gc[ii]);
                    if (tkHITS > temp_hits2) {
                      temp_hits2 = tkHITS;
                      for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
                        EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
                        if (!( * itTrk) -> isEmcShowerValid()) {
                          continue;
                        }
                        RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
                        if (emcTrk -> getShowerId().get_value() == tracks2[z] -> getShowerId().get_value()) {
                          MCtkID2 = ( * itTrk) -> trackId();
                        }
                      }
                    }
                  }
                  if (MCtkID2 == tkID2) {
                    ismatched = true;
                  }
                } // matched the 2nd daughter
              } // Loop over the 2nd daughter
            }
          }
        }
      }
    }
  }

  return ismatched;
}
bool DStarID::ETA3PI(int tkID1, int tkID2, int shID1, int shID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  if (!Par2SH(111, 22, 22, shID1, shID2)) {
    return ismatched;
  }

  int eta_tk1 = -2;
  int eta_tk2 = -2;
  int eta_sh1 = -2;
  int eta_sh2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  int MCID_sh1 = -2;
  int MCID_sh2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 221 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          eta_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 221 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          eta_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 111 &&
          particles[i] -> mother().mother().particleProperty() == 221 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          eta_sh1 = particles[i] -> mother().mother().trackIndex();
          MCID_sh1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID2) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_sh1 &&
          particles[i] -> mother().particleProperty() == 111 &&
          particles[i] -> mother().mother().particleProperty() == 221 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          eta_sh2 = particles[i] -> mother().mother().trackIndex();
          MCID_sh2 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (eta_tk1 != -2 &&
    eta_tk1 == eta_tk2 && eta_tk1 == eta_sh1 && eta_tk1 == eta_sh2 &&
    eta_tk2 == eta_sh1 && eta_tk2 == eta_sh2 &&
    eta_sh1 == eta_sh2) {
    ismatched = true;
  }

  return ismatched;
}
bool DStarID::ETAPGG(int tkID1, int tkID2, int shID1, int shID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  if (!Par2SH(221, 22, 22, shID1, shID2)) {
    return ismatched;
  }

  int etap_tk1 = -2;
  int etap_tk2 = -2;
  int etap_sh1 = -2;
  int etap_sh2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  int MCID_sh1 = -2;
  int MCID_sh2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 331 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 331 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 221 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh1 = particles[i] -> mother().mother().trackIndex();
          MCID_sh1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID2) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_sh1 &&
          particles[i] -> mother().particleProperty() == 221 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh2 = particles[i] -> mother().mother().trackIndex();
          MCID_sh2 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (etap_tk1 != -2 &&
    etap_tk1 == etap_tk2 && etap_tk1 == etap_sh1 && etap_tk1 == etap_sh2 &&
    etap_tk2 == etap_sh1 && etap_tk2 == etap_sh2 &&
    etap_sh1 == etap_sh2) {
    ismatched = true;
  }

  return ismatched;
}
bool DStarID::ETAP3P(int tkID1, int tkID2, int tkID3, int tkID4, int shID1, int shID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  if (!ETA3PI(tkID3, tkID4, shID1, shID2)) {
    return ismatched;
  }

  int etap_tk1 = -2;
  int etap_tk2 = -2;
  int etap_tk3 = -2;
  int etap_tk4 = -2;
  int etap_sh1 = -2;
  int etap_sh2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  int MCID_tk3 = -2;
  int MCID_tk4 = -2;
  int MCID_sh1 = -2;
  int MCID_sh2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 331 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 331 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID3) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> trackIndex() != MCID_tk2 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> mother().particleProperty() == 221 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk3 = particles[i] -> mother().mother().trackIndex();
          MCID_tk3 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID4) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> trackIndex() != MCID_tk2 &&
          particles[i] -> trackIndex() != MCID_tk3 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> mother().particleProperty() == 221 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk4 = particles[i] -> mother().mother().trackIndex();
          MCID_tk4 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 111 &&
          particles[i] -> mother().mother().particleProperty() == 221 &&
          particles[i] -> mother().mother().mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh1 = particles[i] -> mother().mother().mother().trackIndex();
          MCID_sh1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID2) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_sh1 &&
          particles[i] -> mother().particleProperty() == 111 &&
          particles[i] -> mother().mother().particleProperty() == 221 &&
          particles[i] -> mother().mother().mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh2 = particles[i] -> mother().mother().mother().trackIndex();
          MCID_sh2 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (etap_tk1 != -2 &&
    etap_tk1 == etap_tk2 && etap_tk1 == etap_tk3 && etap_tk1 == etap_tk4 && etap_tk1 == etap_sh1 && etap_tk1 == etap_sh2 &&
    etap_tk2 == etap_tk3 && etap_tk2 == etap_tk4 && etap_tk2 == etap_sh1 && etap_tk2 == etap_sh2 &&
    etap_tk3 == etap_tk4 && etap_tk3 == etap_sh1 && etap_tk3 == etap_sh2 &&
    etap_tk4 == etap_sh1 && etap_tk4 == etap_sh2 &&
    etap_sh1 == etap_sh2) {
    ismatched = true;
  }

  return ismatched;
}
bool DStarID::RHOPP(int tkID1, int tkID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;

  int rho_tk1 = -2;
  int rho_tk2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 113 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          rho_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 113 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          rho_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (rho_tk1 != -2 && rho_tk1 == rho_tk2) {
    ismatched = true;
  }

  return ismatched;
}
bool DStarID::ETAPRG(int tkID1, int tkID2, int shID1) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;
  if (!RHOPP(tkID1, tkID2)) {
    return ismatched;
  }

  int etap_tk1 = -2;
  int etap_tk2 = -2;
  int etap_sh1 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  int MCID_sh1 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> mother().particleProperty() == 113 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk1 = particles[i] -> mother().mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().mother().particleProperty() == 331 &&
          particles[i] -> mother().particleProperty() == 113 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          etap_tk2 = particles[i] -> mother().mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == shID1) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 331 &&
          particles[i] -> particleProperty() == 22) {
          temp_hits = relevance;
          etap_sh1 = particles[i] -> mother().trackIndex();
          MCID_sh1 = particles[i] -> trackIndex();
        }
      }
    }
  }

  if (etap_tk1 != -2 &&
    etap_tk1 == etap_tk2 && etap_tk1 == etap_sh1 &&
    etap_tk2 == etap_sh1) {
    ismatched = true;
  }

  return ismatched;
}
bool DStarID::KSPP(int tkID1, int tkID2) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  bool ismatched = false;

  int ks_tk1 = -2;
  int ks_tk2 = -2;
  int MCID_tk1 = -2;
  int MCID_tk2 = -2;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID1) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> mother().particleProperty() == 310 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          ks_tk1 = particles[i] -> mother().trackIndex();
          MCID_tk1 = particles[i] -> trackIndex();
        }
      }
    }
  }
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isMdcKalTrackValid()) {
      continue;
    }
    RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
    if (( * itTrk) -> trackId() == tkID2) {
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits &&
          particles[i] -> trackIndex() != MCID_tk1 &&
          particles[i] -> mother().particleProperty() == 310 &&
          fabs(particles[i] -> particleProperty()) == 211) {
          temp_hits = relevance;
          ks_tk2 = particles[i] -> mother().trackIndex();
          MCID_tk2 = particles[i] -> trackIndex();
        }
      }
    }
  }
  if (ks_tk1 != -2 && ks_tk1 == ks_tk2) {
    ismatched = true;
  }

  return ismatched;
}

bool DStarID::MatchedTK(int tkID, int mcPDG) {
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  bool ismatched = false;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        if (!particles[i] -> decayFromGenerator()) {
          continue;
        }
        if (particles[i] -> decayInFlight()) {
          continue;
        }
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
        }
        if (thebestmatchedPDG == mcPDG) {
          ismatched = true;
        }
      }
    }
  }

  return ismatched;
}

int DStarID::MatchedTKID(int mcID, RecMdcKalTrackVector tracks) {
  int tkID = 0;
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  if (!navigator) {
    //cout << "EventNavigator not found" << endreq;
    tkID = -1;
    return int(tkID);
  }

  for (; iter_mc != mcParticles -> end(); iter_mc++) {
    if (mcID == ( * iter_mc) -> trackIndex()) {
      int temp_tkHITS = 0;
      for (unsigned int z = 0; z < tracks.size(); z++) {
        int tkHITS = navigator -> getMcParticleRelevance(tracks[z], * iter_mc);
        if (tkHITS > temp_tkHITS) {
          temp_tkHITS = tkHITS;
          tkID = tracks[z] -> trackId();
        }
      } // matched tracks
    }
  } // iter_mc

  return int(tkID);
}
int DStarID::MatchedSHID(int mcID, RecEmcShowerVector showers) {
  int shID = 0;
  SmartDataPtr < EvtRecEvent > evtRecEvent(eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < Event::McParticleCol > mcParticles(eventSvc(), EventModel::MC::McParticleCol);
  Event::McParticleCol::iterator iter_mc = mcParticles -> begin();
  SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
  if (!navigator) {
    //cout << "EventNavigator not found" << endreq;
    shID = -1;
    return int(shID);
  }
  int shIDtemp = -1;
  for (; iter_mc != mcParticles -> end(); iter_mc++) {
    if (mcID == ( * iter_mc) -> trackIndex()) {
      int temp_shHITS = 0;
      for (unsigned int z = 0; z < showers.size(); z++) {
        int shHITS = navigator -> getMcParticleRelevance(showers[z], * iter_mc);
        if (shHITS > temp_shHITS) {
          temp_shHITS = shHITS;
          shIDtemp = showers[z] -> getShowerId().get_value();
        }
      } // matched showers
    }
  } // iter_mc

  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (emcTrk -> getShowerId().get_value() == shIDtemp) {
      shID = ( * itTrk) -> trackId();
    }
  }

  return int(shID);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
std::string DStarID::ParticleName(int pdg) {
  IPartPropSvc * p_PartPropSvc;
  static
  const bool CREATEIFNOTTHERE(true);
  StatusCode PartPropStatus = service("PartPropSvc", p_PartPropSvc, CREATEIFNOTTHERE);
  if (!PartPropStatus.isSuccess() || 0 == p_PartPropSvc) {
    std::cout << " Could not initialize Particle Properties Service" << std::endl;
    return "0";
  }
  m_particleTable = p_PartPropSvc -> PDT();

  std::string name;
  if (m_particleTable -> particle(pdg))
    name = m_particleTable -> particle(pdg) -> name();
  else if (m_particleTable -> particle(-pdg))
    name = "anti " + m_particleTable -> particle(-pdg) -> name();

  return name;
}
void DStarID::DumpTree() {

  SmartDataPtr < Event::McParticleCol > mcParticleCol(eventSvc(),
    EventModel::MC::McParticleCol);
  if (mcParticleCol) {
    cout << "-------------------------" << endl;
    //cout << "Run- " << m_run << ", Event- " << m_event << endl;

    //////////////////////////////
    /// Dump vertices
    cout << "Vertices- " << endl;

    // "Cluster" is first particle before BesEvtGen particle
    bool foundClusterAsMother = false;
    //if(!m_BesEvtGenOnly) foundClusterAsMother = true;
    Event::McParticleCol::iterator iter_mc = mcParticleCol -> begin();
    for (; iter_mc != mcParticleCol -> end(); iter_mc++) {

      if (!( * iter_mc) -> primaryParticle()) {
        if (( * iter_mc) -> mother().particleProperty() == 91) foundClusterAsMother = true;
      }
      // took this out for qqbar MC
      //if(!foundClusterAsMother) continue;

      const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();

      if (gc.size() > 0) {
        cout << " " << ParticleName(( * iter_mc) -> particleProperty()) << " [" <<
          ( * iter_mc) -> trackIndex() << "] -> ";

        for (unsigned int ii = 0; ii < gc.size(); ii++) {
          if (ii != (gc.size() - 1))
            cout << ParticleName(gc[ii] -> particleProperty()) << " [" <<
            gc[ii] -> trackIndex() << "], ";
          else
            cout << ParticleName(gc[ii] -> particleProperty()) <<
            " [" << gc[ii] -> trackIndex() << "]" << endl;
        }

      } // End of "gc.size() > 0" IF

    } // End of "McParticleCol" FOR LOOP

    //////////////////////////////////////
    /// Dump particles

    cout << endl << "Particles-  [#Children, primParticle, leafParticle, decayFromGen, decayInFlight] " <<
      endl;

    foundClusterAsMother = false;
    //if(!m_BesEvtGenOnly) foundClusterAsMother = true;

    bool firstDecayInFlight = true;
    iter_mc = mcParticleCol -> begin();
    for (; iter_mc != mcParticleCol -> end(); iter_mc++) {

      if (!( * iter_mc) -> primaryParticle()) {
        if (( * iter_mc) -> mother().particleProperty() == 91) foundClusterAsMother = true;
      }
      // took this out for qqbar MC
      //if(!foundClusterAsMother) continue;

      const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
      int numChildren = gc.size();

      string primaryParticle = "F";
      if (( * iter_mc) -> primaryParticle()) primaryParticle = "T";

      string leafParticle = "F";
      if (( * iter_mc) -> leafParticle()) leafParticle = "T";

      string decayFromGen = "F";
      if (( * iter_mc) -> decayFromGenerator()) decayFromGen = "T";

      string decayInFlight = "F";
      if (( * iter_mc) -> decayInFlight()) {
        decayInFlight = "T";

        if (firstDecayInFlight) {
          cout << endl;
          firstDecayInFlight = false;
        }
      }

      //if(!(*iter_mc)->decayFromGenerator()) 
      //{
      cout << " " << ( * iter_mc) -> trackIndex() << "- " <<
        ParticleName(( * iter_mc) -> particleProperty()) <<
        "  ID = " << ( * iter_mc) -> particleProperty() <<
        " p4 = " << ( * iter_mc) -> initialFourMomentum() << " [" <<
        numChildren << ", " <<
        primaryParticle << ", " << leafParticle << ", " <<
        decayFromGen << ", " << decayInFlight << "]" <<
        endl;
      //}
    } // End of "McParticleCol" FOR LOOP

  }
  cout << endl << endl;
}
int DStarID::countPar(int parentID, int targetID) {
  SmartDataPtr < Event::McParticleCol > mcParticleCol(eventSvc(),
    EventModel::MC::McParticleCol);
  int numcount = 0;
  Event::McParticleCol::iterator iter_mc = mcParticleCol -> begin();
  for (; iter_mc != mcParticleCol -> end(); iter_mc++) {
    if (( * iter_mc) -> particleProperty() == parentID &&
      ( * iter_mc) -> decayFromGenerator() && !( * iter_mc) -> decayInFlight()) {
      const SmartRefVector < Event::McParticle > & gc = ( * iter_mc) -> daughterList();
      if (gc.size() > 0) {
        for (unsigned int ii = 0; ii < gc.size(); ii++) {
          if (gc[ii] -> particleProperty() == targetID &&
            gc[ii] -> decayFromGenerator() && !gc[ii] -> decayInFlight()) {
            numcount++;
          }
          if (gc[ii] -> daughterList().size() > 0) {
            const SmartRefVector < Event::McParticle > & ggc = gc[ii] -> daughterList();
            for (unsigned int iii = 0; iii < ggc.size(); iii++) {
              if (ggc[iii] -> particleProperty() == targetID &&
                ggc[iii] -> decayFromGenerator() && !ggc[iii] -> decayInFlight()) {
                numcount++;
              }
              if (ggc[iii] -> daughterList().size() > 0) {
                const SmartRefVector < Event::McParticle > & gggc = ggc[iii] -> daughterList();
                for (unsigned int iiii = 0; iiii < gggc.size(); iiii++) {
                  if (gggc[iiii] -> particleProperty() == targetID &&
                    gggc[iiii] -> decayFromGenerator() && !gggc[iiii] -> decayInFlight()) {
                    numcount++;
                  }
                  if (gggc[iiii] -> daughterList().size() > 0) {
                    const SmartRefVector < Event::McParticle > & ggggc = gggc[iiii] -> daughterList();
                    for (unsigned int iiiii = 0; iiiii < ggggc.size(); iiiii++) {
                      if (ggggc[iiiii] -> particleProperty() == targetID &&
                        ggggc[iiiii] -> decayFromGenerator() && !ggggc[iiiii] -> decayInFlight()) {
                        numcount++;
                      }
                      if (ggggc[iiiii] -> daughterList().size() > 0) {
                        const SmartRefVector < Event::McParticle > & gggggc = ggggc[iiiii] -> daughterList();
                        for (unsigned int iiiiii = 0; iiiiii < gggggc.size(); iiiiii++) {
                          if (gggggc[iiiiii] -> particleProperty() == targetID &&
                            gggggc[iiiiii] -> decayFromGenerator() && !gggggc[iiiiii] -> decayInFlight()) {
                            numcount++;
                          }
                          if (gggggc[iiiiii] -> daughterList().size() > 0) {
                            const SmartRefVector < Event::McParticle > & ggggggc = gggggc[iiiiii] -> daughterList();
                            for (unsigned int iiiiiii = 0; iiiiiii < ggggggc.size(); iiiiiii++) {
                              if (ggggggc[iiiiiii] -> particleProperty() == targetID &&
                                ggggggc[iiiiiii] -> decayFromGenerator() && !ggggggc[iiiiiii] -> decayInFlight()) {
                                numcount++;
                              }
                            } // iiiiiii
                          }
                        } // iiiiii
                      }
                    } // iiiii
                  }
                } // iiii
              }
            } // iii
          }
        } // ii
      }
    } // parent found
  }
  return numcount;
}
bool DStarID::isTagTrue(DTagToolIterator iter_dtag) {
    bool isfound = false;
    SmartDataPtr < EventNavigator > navigator(eventSvc(), "/Event/Navigator");
    if (!navigator) {
      return isfound;
    }
    int m_mode = ( * iter_dtag) -> decayMode();
    int m_charge = ( * iter_dtag) -> charge();
    SmartRefVector < EvtRecTrack > tracks = ( * iter_dtag) -> tracks();
    SmartRefVector < EvtRecTrack > showers = ( * iter_dtag) -> showers();
    vector < int > numchan;
    numchan.clear();
    //K+/- :1 
    //Pi+/-:2
    //Pi0  :3
    //Eta  :4
    //Ks   :5
    //eta'(pipieta)         :6
    //eta'(rhogamma)        :7
    //eta'(pipieta(pipipi0)):8
    //eta(pipipi0)          :9

    //cout << " *********** the mode = " << m_mode << endl;

    //if (m_mode != 404) {return isfound;}

    if (m_mode == 400) {
      numchan.push_back(5);
      numchan.push_back(1);
    } // DstoKsK
    else if (m_mode == 401) {
      numchan.push_back(1);
      numchan.push_back(1);
      numchan.push_back(2);
    } // DstoKDStar
    else if (m_mode == 402) {
      numchan.push_back(5);
      numchan.push_back(1);
      numchan.push_back(3);
    } // DstoKsDStar0
    else if (m_mode == 403) {
      numchan.push_back(5);
      numchan.push_back(5);
      numchan.push_back(2);
    } // DstoKsKsPi
    else if (m_mode == 404) {
      numchan.push_back(1);
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(3);;
    } // DstoKDStarPi0
    else if (m_mode == 405) {
      numchan.push_back(5);
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoKsKplusPiPi
    else if (m_mode == 406) {
      numchan.push_back(5);
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoKsKminusPiPi
    else if (m_mode == 407) {
      numchan.push_back(1);
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoKDStarPiPi
    else if (m_mode == 420) {
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoPiPi0
    else if (m_mode == 421) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoPiPiPi
    else if (m_mode == 422) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoPiPiPiPi0
    else if (m_mode == 423) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoPiPiPiPiPi
    else if (m_mode == 424) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoPiPiPiPiPiPi0
    else if (m_mode == 425) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(3);
    } // DstoPiPiPiPi0Pi0 <-- this is wrong
    else if (m_mode == 440) {
      numchan.push_back(2);
      numchan.push_back(4);
    } // DstoPiEta
    else if (m_mode == 441) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(4);
    } // DstoPiPi0Eta
    else if (m_mode == 442) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(4);
    } // DstoPiPiPiEta
    else if (m_mode == 450) {
      numchan.push_back(2);
      numchan.push_back(9);
    } // DstoPiEtaPiPiPi0
    else if (m_mode == 451) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(9);
    } // DstoPiPi0EtaPiPiPi0
    else if (m_mode == 452) {
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(9);
    } // DstoPiPiPiEtaPiPiPi0
    else if (m_mode == 460) {
      numchan.push_back(2);
      numchan.push_back(6);
    } // DstoPiEPPiPiEta
    else if (m_mode == 461) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(6);
    } // DstoPiPi0EPPiPiEta
    else if (m_mode == 470) {
      numchan.push_back(2);
      numchan.push_back(8);
    } // DstoPiEPPiPiEtaPiPiPi0
    else if (m_mode == 471) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(8);
    } // DstoPiPi0EPPiPiEtaPiPiPi0
    else if (m_mode == 480) {
      numchan.push_back(2);
      numchan.push_back(7);
    } // DstoPiEPRhoGam
    else if (m_mode == 481) {
      numchan.push_back(2);
      numchan.push_back(3);
      numchan.push_back(7);
    } // DstoPiPi0EPRhoGam
    else if (m_mode == 500) {
      numchan.push_back(5);
      numchan.push_back(2);
    } // DstoKsPi
    else if (m_mode == 501) {
      numchan.push_back(5);
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoKsPiPi0
    else if (m_mode == 502) {
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
    } // DstoDStarPi
    else if (m_mode == 503) {
      numchan.push_back(1);
      numchan.push_back(2);
      numchan.push_back(2);
      numchan.push_back(3);
    } // DstoDStarPiPi0
    else if (m_mode == 504) {
      numchan.push_back(1);
      numchan.push_back(1);
      numchan.push_back(1);
    } // DstoKKK
    bool isparmatched = true;
    int npi = 0;
    int nka = 0;
    int np0 = 0;
    int net = 0;
    int nks = 0;
    int nep6 = 0;
    int nep7 = 0;
    int nep8 = 0;
    int net9 = 0;
    for (int i = 0; i < ( * iter_dtag) -> numOfChildren(); i++) {
      //cout << " EACH PARTICLES = " << numchan[i] << endl;
      if (numchan[i] == 1) { // Kaon
        nka++;
        if (!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 321, 0, 431)) &&
          !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 321, 0, -431))) {
          isparmatched = false;
        }
        //cout << " K1 = " << MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),321,0, 431) << endl
        //   << " K2 = " << MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),321,0,-431) << endl
        //   << " K matched = " << isparmatched << endl;
      } else if (numchan[i] == 2) { // pion
        npi++;
        if (!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 0, 431)) &&
          !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 0, -431))) {
          isparmatched = false;
        }
        //cout << " P1 = " << MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,0, 431) << endl
        //   << " P2 = " << MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,0,-431) << endl
        //   << " P matched = " << isparmatched << endl;
      } else if (numchan[i] == 3) { // pi0->gg
        np0++;
        //if ( (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431)) && 
        //    !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431))) ||
        //   (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111, 431)) && 
        //    !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111,-431))) ) {isparmatched=false;}

        //cout << " P01 = " << (m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431)) << endl
        //   << " P02 = " << (m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431)) << endl
        //   << " P03 = " << (m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111, 431)) << endl
        //   << " P04 = " << (m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111,-431)) << endl
        //   << " P0 matched = " << isparmatched << endl;
        //cout << " Try1 = " << (m_charge>0 &&
        //		       (MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431) ||
        //			MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111, 431))) << endl
        //   << " Try2 = " << (m_charge<0 &&
        //		       (MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431) ||
        //			MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111,-431))) << endl;

        if (!(m_charge > 0 &&
            (MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 111, 431) ||
              MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 111, 431))) &&
          !(m_charge < 0 &&
            (MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 111, -431) ||
              MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 111, -431)))) {
          isparmatched = false;
        }

      } else if (numchan[i] == 4) { // eta->gg
        net++;
        //if ( (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,221, 431)) && 
        //    !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,221,-431))) ||
        //   (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,221, 431)) && 
        //    !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,221,-431))) ) {isparmatched=false;} 

        if (!(m_charge > 0 &&
            (MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 221, 431) ||
              MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 221, 431))) &&
          !(m_charge < 0 &&
            (MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 221, -431) ||
              MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 221, -431)))) {
          isparmatched = false;
        }
      } else if (numchan[i] == 5) { // Ks
        nks++;
        if ((!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 2] -> trackId(), 211, 310, 431)) &&
            !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 2] -> trackId(), 211, 310, -431))) ||
          (!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 310, 431)) &&
            !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 310, -431)))) {
          isparmatched = false;
        }
      } else if (numchan[i] == 6) { // eta' -> pipi eta(->gg)
        nep6++;
        if ((!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 2] -> trackId(), 211, 331, 431)) &&
            !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 2] -> trackId(), 211, 331, -431))) ||
          (!(m_charge > 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 331, 431)) &&
            !(m_charge < 0 && MCTKPID(tracks[nka + npi + 2 * nks + 2 * nep6 + 2 * nep7 + 4 * nep8 + 2 * net9 - 1] -> trackId(), 211, 331, -431))) ||
          (!(m_charge > 0 && MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 221, 431)) &&
            !(m_charge < 0 && MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 2] -> trackId(), 22, 221, -431))) ||
          (!(m_charge > 0 && MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 221, 431)) &&
            !(m_charge < 0 && MCSHPID(showers[2 * np0 + 2 * net + 2 * nep6 + nep7 + 2 * nep8 + 2 * net9 - 1] -> trackId(), 22, 221, -431)))) {
          isparmatched = false;
        }
      } else if (numchan[i] == 7) { // eta' -> rho gamma
                nep7++;
        if ( (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,113, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,113,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,113, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,113,-431))) ||
             (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,331, 431)) &&
              !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,331,-431))) ) {isparmatched=false;} }
      else if(numchan[i]==8){ // eta' -> pipi eta(->pipipi0)
        nep8++;

        if ( (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-4]->trackId(),211,331, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-4]->trackId(),211,331,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-3]->trackId(),211,331, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-3]->trackId(),211,331,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,221, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,221,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,221, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,221,-431))) ||
             (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431)) &&
              !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431))) ) {isparmatched=false;} }
      else if(numchan[i]==9){ // eta -> pipipi0
        net9++;
        if ( (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,221, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-2]->trackId(),211,221,-431))) ||
             (!(m_charge>0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,221, 431)) &&
              !(m_charge<0&&MCTKPID(tracks[nka+npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9-1]->trackId(),211,221,-431))) ||
             (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111, 431)) &&
              !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-2]->trackId(),22,111,-431))) ||
             (!(m_charge>0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111, 431)) &&
              !(m_charge<0&&MCSHPID(showers[2*np0+2*net+2*nep6+nep7+2*nep8+2*net9-1]->trackId(),22,111,-431))) ) {isparmatched=false;} }
      //if (m_mode == 404)
      // {
      //   cout << " EACH PARTICLES = " << numchan[i] << endl;
      // }
    }
      if (!isparmatched) {return isfound;}

  //cout << "============ B "<< isparmatched << endl;

  int totalPi = npi+2*nks+2*nep6+2*nep7+4*nep8+2*net9;
  int totalKa = nka;
  int totalP0 = np0+nep8+net9;
  int totalET = net+nep6+nep8+net9;
  int totalKs = nks;
  int totalRo = nep7;
  int totalEP = nep6+nep7+nep8;
  if ( m_charge>0 &&
       (countPar( 431, 211)+countPar( 431,-211))==totalPi &&
       (countPar( 431, 321)+countPar( 431,-321))==totalKa &&
       countPar( 431,111)==totalP0&&countPar( 431, 221)==totalET&&
       countPar( 431,310)==totalKs&&countPar( 431, 113)==totalRo&&countPar( 431,331)==totalEP ) {isfound=true;}
  if ( m_charge<0 &&
       (countPar(-431, 211)+countPar(-431,-211))==totalPi &&
       (countPar(-431, 321)+countPar(-431,-321))==totalKa &&
       countPar(-431,111)==totalP0&&countPar(-431, 221)==totalET&&
       countPar(-431,310)==totalKs&&countPar(-431, 113)==totalRo&&countPar(-431,331)==totalEP ) {isfound=true;}

  //cout << "np0 = " << np0 << " nep8 = " << nep8 << " net0 = " << net9 << endl;
  //cout << "============ C "<< isfound << endl;

  //if(m_mode==450 && m_charge>0)
  //if(m_mode==401)
  //if (m_charge>0)
  //{
  //  cout << "------------------------" << endl;
  //  cout << " mode = " << m_mode << endl;
  //  cout << " pi+-= " << countPar( 431, 211)+countPar( 431,-211) << " expected = " << totalPi << endl
  //       << " K+- = " << countPar( 431, 321)+countPar( 431,-321) << " expected = " << totalKa << endl
  //       << " pi0 = " << countPar( 431, 111) << " expected = " << totalP0 << endl
  //       << " eta = " << countPar( 431, 221) << " expected = " << totalET << endl
  //       << " Ks  = " << countPar( 431, 310) << " expected = " << totalKs << endl
  //       << " Ro  = " << countPar( 431, 113) << " expected = " << totalRo << endl
  //       << " EP  = " << countPar( 431, 331) << " expected = " << totalEP << endl;
  //}

  return isfound;
}
