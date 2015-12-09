'''

Run LFV H->MuTau analysis in the mu+tau channel.

Authors: Maria Cepeda, Aaron Levine, Evan K. Friis, UW

'''

import MuMuTauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
#import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
#import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import ROOT

#requireZeroJets = bool ('true' in os.environ['zeroJet']) #is NUP==5
#print "env zeroJet:", os.environ['zeroJet']
#print "require zeroJets:", requireZeroJets

#isEmbed = bool ('true' in os.environ['embed'])
#print "isEmbed:", isEmbed

################################################################################
#### FUNCTION DEFINITION #######################################################
################################################################################

from math import sqrt, pi

def deltaPhi(phi1, phi2):
  PHI = abs(phi1-phi2)
  if PHI<=pi:
      return PHI
  else:
      return 2*pi-PHI


################################################################################
#### MC-DATA and PU corrections ################################################
################################################################################

# Determine MC-DATA corrections
#is7TeV = bool('7TeV' in os.environ['jobid'])
#print "Is 7TeV:", is7TeV

# Make PU corrector from expected data PU distribution
# PU corrections .root files from pileupCalc.py
#pu_distributions = glob.glob(os.path.join(
#    'inputs', os.environ['jobid'], 'data_TauPlusX*pu.root'))
#	'inputs', os.environ['jobid'], 'data_DoubleMu*2012*pu.root'))
#pu_corrector = PileupWeight.PileupWeight(
#    'S6' if is7TeV else 'S10', *pu_distributions)

#muon_pog_PFTight_2012 = MuonPOGCorrections.make_muon_pog_PFTight_2012()
#muon_pog_PFRelIsoDB02_2012 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB012_2012()
#muon_pog_IsoMu24eta2p1_2012 = MuonPOGCorrections.make_muon_pog_IsoMu24eta2p1_2012()


# Get object ID and trigger corrector functions
def mc_corrector_2012(row):
    if row.run > 2:
        return 1
#    pu = pu_corrector(row.nTruePU)
#    m1id = muon_pog_PFTight_2012(row.m1Pt, row.m1Eta)
#    m1iso = muon_pog_PFRelIsoDB02_2012(row.m1Pt, row.m1Eta)
#    m_trg = H2TauCorrections.correct_mueg_mu_2012(row.m1Pt, row.m1AbsEta)
#    m_trg = muon_pog_IsoMu24eta2p1_2012(row.m1Pt, row.m1AbsEta)
#    return pu*m1id*m1iso*m_trg
    return 1

# Determine which set of corrections to use
mc_corrector = mc_corrector_2012

################################################################################
#### HISTOGRAM BOOKING ######## ################################################
################################################################################

class AnalyzeMuMuTauTight_FakeRate(MegaBase):
    tree = 'mmt/final/Ntuple'

    def __init__(self, tree, outfile, **kwargs):
        super(AnalyzeMuMuTauTight_FakeRate, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuMuTauTree.MuMuTauTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):
        self.book("Yield", "looseVsPt_old", "", 200, 0, 200)
        self.book("Yield", "looseVsEta_old", "", 200, -2.5, 2.5)
        self.book("Yield", "looseVstDM_old", "" , 11,  -0.5 , 10.5)

        self.book("Yield", "tightVsPt_old", "", 200, 0, 200)
        self.book("Yield", "tightVsEta_old", "", 200, -2.5, 2.5)
        self.book("Yield", "tightVstDM_old", "" , 11,  -0.5 , 10.5)

        self.book("Yield", "vlooseVsPt_new", "", 200, 0, 200)
        self.book("Yield", "vlooseVsEta_new", "", 200, -2.5, 2.5)
        self.book("Yield", "vlooseVstDM_new", "" , 11,  -0.5 , 10.5)

        self.book("Yield", "looseVsPt_new", "", 200, 0, 200)
        self.book("Yield", "looseVsEta_new", "", 200, -2.5, 2.5)
        self.book("Yield", "looseVstDM_new", "" , 11,  -0.5 , 10.5)

        self.book("Yield", "mediumVsPt_new", "", 200, 0, 200)
        self.book("Yield", "mediumVsEta_new", "", 200, -2.5, 2.5)
        self.book("Yield", "mediumVstDM_new", "" , 11,  -0.5 , 10.5)

        self.book("Yield", "tightVsPt_new", "", 200, 0, 200)
        self.book("Yield", "tightVsEta_new", "", 200, -2.5, 2.5)
        self.book("Yield", "tightVstDM_new", "" , 11,  -0.5 , 10.5)
        
        self.book("Yield", "vtightVsPt_new", "", 200, 0, 200)
        self.book("Yield", "vtightVsEta_new", "", 200, -2.5, 2.5)
        self.book("Yield", "vtightVstDM_new", "" , 11,  -0.5 , 10.5)

        self.book("Yield", "vvtightVsPt_new", "", 200, 0, 200)
        self.book("Yield", "vvtightVsEta_new", "", 200, -2.5, 2.5)
        self.book("Yield", "vvtightVstDM_new", "" , 11,  -0.5 , 10.5)

    def correction(self, row):
        return mc_corrector(row)

    def fill_yield(self, row, histo, value):
        histos = self.histograms
        weight = self.correction(row)
        histos['Yield/'+histo].Fill(value, weight)

################################################################################
#### SELECTIONS ################################################################
################################################################################
        
    def presel(self, row):
      #if not row.mu17mu8Pass:
      #  return False
      if row.m1_m2_Mass < 80 or row.m1_m2_Mass > 100:
        return False
      return True

    def presel_vbf(self, row):
#      if row.jetVeto30 != 2:
#        return False
      return True

    def kinematics(self, row):
        if row.m1Pt < 25:
            return False
        if row.m2Pt < 25:
            return False
        if row.tPt < 20:
            return False
        if abs(row.m1Eta) >= 2.1:
            return False
        if abs(row.m2Eta) >= 2.1:
            return False
        if abs(row.tEta)>=2.3 :
            return False

        return True

    def gg0(self,row):
      if row.m1Pt < 45:
        return False
      if row.tPt < 35:
        return False
      
      if row.tMtToPfMet_Ty1>50:
        return False
#      if row.mMtToPfMet_Ty1<30:
#        return False
      
#      if (deltaPhi(row.tPhi,row.pfMetPhi)>0.3):
#        return False
#      if (deltaPhi(row.mPhi,row.pfMetPhi)<2.5):
#        return False        
      if (deltaPhi(row.m1Phi,row.tPhi)<2.7):    
        return False
      return True    

    def vbf(self,row):
      if abs(row.vbfDeta)<3.5:
        return False
      if row.vbfMass<550:
        return False
      if row.vbfJetVeto30>0:
        return False
      if row.tMtToPfMet_Ty1>35:
        return False
      return True

    def oppositesign(self,row):
      if row.m1Charge*row.tCharge!=-1:
        return False
      return True

    def samesign(self, row):
      if row.m1Charge != row.tCharge:
        return False
      return True

    def tDecayMode(self, row):
      return row.tDecayMode

    def obj1_id(self, row):
      return bool(row.m1PFIDTight)  and bool(abs(row.m1PVDZ) < 0.2) and bool(row.m2PFIDTight)  and bool(abs(row.m2PVDZ) < 0.2) 

    def obj2_idOld(self, row):
      return  row.tAgainstElectronLooseMVA5 and row.tAgainstMuonTight3 and row.tDecayModeFinding
                  
    def obj2_idNew(self, row):
      return  row.tAgainstElectronLooseMVA5 and row.tAgainstMuonTight3 and row.tDecayModeFindingNewDMs
    
    def vetos(self,row):
      return bool (row.muVetoPt5IsoIdVtx<1) and bool (row.tauVetoPt20Loose3HitsVtx<1)
    #bool (row.eVetoCicTightIso<1) and 
    #bool (row.bjetVeto<1)

    def obj1_iso(self, row):
      return bool(row.m1RelPFIsoDBDefault <0.12) and bool(row.m2RelPFIsoDBDefault <0.12)

    def obj2_looseIsoOld(self, row):
      return row.tByLooseCombinedIsolationDeltaBetaCorr3Hits
    
    def obj2_tightIsoOld(self, row):
      return row.tByTightCombinedIsolationDeltaBetaCorr3Hits
    
    def obj2_VVTightIsoNew(self, row):
      return row.tByVVTightIsolationMVA3newDMwLT
    def obj2_VTightIsoNew(self, row):
      return row.tByVTightIsolationMVA3newDMwLT
    def obj2_TightIsoNew(self, row):
      return row.tByTightIsolationMVA3newDMwLT
    def obj2_MediumIsoNew(self, row):
      return row.tByMediumIsolationMVA3newDMwLT
    def obj2_LooseIsoNew(self, row):
      return row.tByLooseIsolationMVA3newDMwLT
    def obj2_VLooseIsoNew(self, row):
      return row.tByVLooseIsolationMVA3newDMwLT


    def obj1_antiiso(self, row):
      return bool(row.m1RelPFIsoDBDefault >0.2) 

    def obj2_antiiso(self, row):
      return  not row.tByTightCombinedIsolationDeltaBetaCorr3Hits

################################################################################
#### RUN!    ###################################################################
################################################################################
    
    def process(self):
      event =0
      sel=False
      for row in self.tree:

        if event!=row.evt:   # This is just to ensure we get the (Mu,Tau) with the highest Pt
          event=row.evt    # In principle the code saves all the MU+Tau posibilities, if an event has several combinations
          sel = False      # it will save them all.
          if sel==True:
            continue 
          
        if not self.presel(row):
          continue

        if not self.kinematics(row):
          continue          
        if not self.obj1_id(row): 
          continue
        #if not self.obj2_id(row):
        #  continue
          
        if not self.vetos(row):
          continue

        sel = True

        #if not self.gg0(row):
        #  continue
        
        if not self.obj1_iso(row):
          continue
        if not self.oppositesign(row):
          continue

        #fake rate with old id/iso
        if self.obj2_idOld(row):
          if self.obj2_looseIsoOld(row):
            self.fill_yield(row, 'looseVsPt_old', row.tJetPt)
            self.fill_yield(row, 'looseVsEta_old', row.tEta)
            self.fill_yield(row, 'looseVstDM_old', row.tDecayMode)
          if self.obj2_tightIsoOld(row):
            self.fill_yield(row, 'tightVsPt_old', row.tJetPt)
            self.fill_yield(row, 'tightVsEta_old', row.tEta)
            self.fill_yield(row, 'tightVstDM_old', row.tDecayMode)

        #fake rate with new id/iso
        if self.obj2_idNew(row):
          if self.obj2_VLooseIsoNew(row):
            self.fill_yield(row, 'vlooseVsPt_new', row.tJetPt)
            self.fill_yield(row, 'vlooseVsEta_new', row.tEta)
            self.fill_yield(row, 'vlooseVstDM_new', row.tDecayMode)
          if self.obj2_LooseIsoNew(row):
            self.fill_yield(row, 'looseVsPt_new', row.tJetPt)
            self.fill_yield(row, 'looseVsEta_new', row.tEta)
            self.fill_yield(row, 'looseVstDM_new', row.tDecayMode)
          if self.obj2_MediumIsoNew(row):
            self.fill_yield(row, 'mediumVsPt_new', row.tJetPt)
            self.fill_yield(row, 'mediumVsEta_new', row.tEta)
            self.fill_yield(row, 'mediumVstDM_new', row.tDecayMode)
          if self.obj2_TightIsoNew(row):
            self.fill_yield(row, 'tightVsPt_new', row.tJetPt)
            self.fill_yield(row, 'tightVsEta_new', row.tEta)
            self.fill_yield(row, 'tightVstDM_new', row.tDecayMode)
          if self.obj2_VTightIsoNew(row):
            self.fill_yield(row, 'vtightVsPt_new', row.tJetPt)
            self.fill_yield(row, 'vtightVsEta_new', row.tEta)
            self.fill_yield(row, 'vtightVstDM_new', row.tDecayMode)
          if self.obj2_VVTightIsoNew(row):
            self.fill_yield(row, 'vvtightVsPt_new', row.tJetPt)
            self.fill_yield(row, 'vvtightVsEta_new', row.tEta)
            self.fill_yield(row, 'vvtightVstDM_new', row.tDecayMode)

            
    def finish(self):
        self.write_histos()

