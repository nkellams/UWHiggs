'''

LFV Higgs->MuTau analysis in the mu-tau channel.

Authors: Aaron Levine, Evan K. Friis, Maria Cepeda, UW

'''

import MuTauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
#import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
#import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import ROOT
import math

#NoZtautau = bool ('true' in os.environ['noZtt'])
#isEmbed = bool ('true' in os.environ['embed'])
#isZeroJet = bool ('true' in os.environ['zeroJet']) #jetVeto30==0 NUP==5
#print "env zeroJet:", os.environ['zeroJet']
#print "isZeroJet:", isZeroJet
#print "NoZtautau:", NoZtautau
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

def collMass_type1(row):
  taupx = row.tPt*math.cos(row.tPhi)
  taupy = row.tPt*math.sin(row.tPhi)
  taupt= row.tPt  
  metpx = row.type1_pfMetEt*math.cos(row.type1_pfMetPhi)
  metpy = row.type1_pfMetEt*math.sin(row.type1_pfMetPhi)
  met = row.type1_pfMetEt
  METproj= abs(metpx*taupx+metpy*taupy)/taupt
  xth=taupt/(taupt+METproj)
  den=math.sqrt(xth)
  mass=row.m_t_Mass/den
  return mass

def collMass(row):
  taupx = row.tPt*math.cos(row.tPhi)
  taupy = row.tPt*math.sin(row.tPhi)
  taupt= row.tPt  
  metpx = row.raw_pfMetEt*math.cos(row.raw_pfMetPhi)
  metpy = row.raw_pfMetEt*math.sin(row.raw_pfMetPhi)
  met = row.raw_pfMetEt
  METproj= abs(metpx*taupx+metpy*taupy)/taupt
  xth=taupt/(taupt+METproj)
  den=math.sqrt(xth)
  mass=row.m_t_Mass/den
  return mass

def fullMT(met,mupt,taupt, metphi, muphi, tauphi):
  mux=mupt*math.cos(muphi)
  muy=mupt*math.sin(muphi)
  metx=met*math.cos(metphi)
  mety=met*math.sin(metphi)
  taux=taupt*math.cos(tauphi)
  tauy=taupt*math.sin(tauphi)
  full_pt=met+mupt+taupt
  full_x=metx+mux+taux
  full_y=mety+muy+tauy
  full_mt_2 = full_pt*full_pt-full_x*full_x-full_y*full_y
  full_mt=0
  if (full_mt_2>0):
    full_mt= math.sqrt(full_mt_2)
  return full_mt
  
def fullPT(met,mupt,taupt, metphi, muphi, tauphi):
  mux=mupt*math.cos(muphi)
  muy=mupt*math.sin(muphi)
  metx=met*math.cos(metphi)
  mety=met*math.sin(metphi)
  taux=taupt*math.cos(tauphi)
  tauy=taupt*math.sin(tauphi)
  full_x=metx+mux+taux
  full_y=mety+muy+tauy
  full_pt_2 = full_x*full_x+full_y*full_y
  full_pt=0
  if (full_pt_2>0):
    full_pt= math.sqrt(full_pt_2)
  return full_pt

################################################################################
#### MC-DATA and PU corrections ################################################
################################################################################

# Determine MC-DATA corrections
#is7TeV = bool('7TeV' in os.environ['jobid'])
#if is7TeV: 
#  print "This code is not setup for 7TeV!:", is7TeV

# Make PU corrector from expected data PU distribution
# PU corrections .root files from pileupCalc.py
pu_distributions = glob.glob(os.path.join('inputs', os.environ['jobid'], 'data_SingleMu*pu.root'))
pu_corrector = PileupWeight.PileupWeight('Asympt25ns', *pu_distributions)

#muon_pog_PFTight_2012 = MuonPOGCorrections.make_muon_pog_PFTight_2012()             
#muon_pog_PFRelIsoDB012_2012 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB012_2012() 
#muon_pog_IsoMu24eta2p1_2012 = MuonPOGCorrections.make_muon_pog_IsoMu24eta2p1_2012()

#muon_pog_PFRelIsoDB02_2012 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB012_2012()


# Note: This should be updated, we are using the prompt-reco corrections

def mc_corrector_2012(row):
  if row.run > 2:
    return 1
  pu = pu_corrector(row.nTruePU)
  #m1id = muon_pog_PFTight_2012(row.mPt, row.mEta)  
  #m1iso = muon_pog_PFRelIsoDB02_2012(row.mPt, row.mEta)
  #m_trg =  muon_pog_IsoMu24eta2p1_2012(row.mPt, row.mAbsEta) 
  #return pu*m1id*m1iso*m_trg
  return pu

mc_corrector = mc_corrector_2012

def getFakeRateFactor(row):
  if row.tDecayMode==0:
    return 0.428
  if row.tDecayMode==1:
    return 0.495
  if row.tDecayMode==10:
    return 0.319


# use Maria's table for now
#  if row.tEta>=1.5:
#    if row.tDecayMode==0:
#      return 0.53
#    if row.tDecayMode==1:
#      return 0.52
#    if row.tDecayMode==10:
#      return 0.32
#  if row.tEta<1.5:
#    if row.tDecayMode==0:
#      return 0.45
#    if row.tDecayMode==1:
#      return 0.44
#    if row.tDecayMode==10:
#      return 0.31
  

################################################################################
#### HISTOGRAM BOOKING ######## ################################################
################################################################################


class AnalyzeGG5(MegaBase):
  tree = 'mt/final/Ntuple'
  
  def __init__(self, tree, outfile, **kwargs):
    super(AnalyzeGG5, self).__init__(tree, outfile, **kwargs)
    # Use the cython wrapper
    self.tree = MuTauTree.MuTauTree(tree)
    self.out = outfile
    # Histograms for each category
    self.histograms = {}
    self.is7TeV = '7TeV' in os.environ['jobid']

  def begin(self):
    names=[
           "gg", "gg_ss", "gg_looseNotTightTau", "gg_ssLooseNotTightTau",
           "pre_gg", "pre_gg_ss", "pre_gg_looseNotTightTau", "pre_gg_ssLooseNotTightTau",
           ]



    for x in range(0, len(names)):
        
      self.book(names[x], "weight", "Event weight", 100, 0, 5)
#      self.book(names[x], "weight_nopu", "Event weight without PU", 100, 0, 5)
#      self.book(names[x], "EmbPtWeight", "embedded weight", 100, 0, 1)
        
#      self.book(names[x], "rho", "Fastjet #rho", 100, 0, 25)
#      self.book(names[x], "nvtx", "Number of vertices", 31, -0.5, 30.5)
#      self.book(names[x], "prescale", "HLT prescale", 21, -0.5, 20.5)
        
#      self.book(names[x], "NUP", "LHE NUP", 11, -0.5, 10.5)
#      self.book(names[x], "ZTauTau", "Is ZTauTau?", 2, -0.5, 1.5)
#      self.book(names[x], 'Ztautau_before0', 'is Ztautua event before req 0', 2,-0.5,1.5)
        
      self.book(names[x], "mPt", "Muon  Pt", 500,0,500)
      self.book(names[x], "mEta", "Muon  eta", 200, -2.5, 2.5)
#      self.book(names[x], "mMtToMVAMET", "Muon MT (MVA)", 500,0,500)
#      self.book(names[x], "mMtToPfMet_Ty1", "Muon MT (PF Ty1)", 500,0,500)
#      self.book(names[x], "mCharge", "Muon Charge", 5, -2, 2)
#      self.book(names[x], "mRelPFIsoDBHtt", "Muon Isolation", 1000, 0, 1)
        
      self.book(names[x], "tPt", "Tau  Pt", 500,0,500)
      self.book(names[x], "tJetPt", "Tau  Jet Pt", 500,0,500)
      self.book(names[x], "tEta", "Tau  eta", 200, -2.5, 2.5)
#      self.book(names[x], "tMtToMVAMET", "Tau MT (MVA)", 500,0,500)
      self.book(names[x], "tMtToPfMet_type1", "Tau MT (PF type1)", 500,0,500)
#      self.book(names[x], "tCharge", "Tau  Charge", 5, -2, 2)
      self.book(names[x], "tDecayMode", "Tau Decay Mode", 21, -0.5, 20.5)
#      self.book(names[x], "tLooseIso", "Tau Isolation Loose", 3, 0, 2)
      self.book(names[x], "tByTightCombinedIsolationDeltaBetaCorr3Hits", "Tau Isolation Tight", 3, 0, 2)
      self.book(names[x], "tByLooseCombinedIsolationDeltaBetaCorr3Hits", "Tau Isolation Loose", 3, 0, 2)
        
#      self.book(names[x], 'mPixHits', 'Mu 1 pix hits', 10, -0.5, 9.5)
#      self.book(names[x], 'mJetBtag', 'Mu 1 JetBtag', 100, -5.5, 9.5)
        
#      self.book(names[x], "LT", "ht", 500, 0, 500)
#      self.book(names[x],"fullMT_mva","fullMT_mva",500,0,500); 
#      self.book(names[x],"fullMT_type1","fullMT_type1",500,0,500);
      self.book(names[x],"collMass_type1","collMass_type1",500,0,500);
      self.book(names[x],"collMass","collMass",500,0,500);
      
#      self.book(names[x],"fullPT_mva","fullPT_mva",500,0,500);
#      self.book(names[x],"fullPT_type1","fullPT_type1",500,0,500);

        
      self.book(names[x], "type1_pfMetEt", "Type1 MET", 500,0,500)
      self.book(names[x], "type1_pfMetPhi", "Type1 MET Phi", 100,-3.15,3.15)
#      self.book(names[x], "mva_metEt", "MVA MET", 500,0,500)

      self.book(names[x], "raw_pfMetEt","raw MET", 500,0,500)
      self.book(names[x], "raw_pfMetPhi","raw MET Phi",100,-3.15,3.15)
#      self.book(names[x], "mva_metPhi","MVA MET Phi",100,-3.15,3.15)
        
#      self.book(names[x], "m_t_Mass", "Muon + Tau Mass", 500,0,500)
#      self.book(names[x], "m_t_Pt", "Muon + Tau Pt", 500,0,500)
#      self.book(names[x], "m_t_DR", "Muon + Tau DR", 100, 0, 10)
#      self.book(names[x], "m_t_DPhi", "Muon + Tau DPhi", 100, 0, 4)
#      self.book(names[x], "m_t_SS", "Muon + Tau SS", 5, -2, 2)
#      self.book(names[x], "m_t_ToMETDPhi_type1", "Muon Tau DPhi to MET", 100, 0, 4)
        
      # Vetoes
#      self.book(names[x], 'bjetVeto', 'Number of b-jets', 5, -0.5, 4.5)
#      self.book(names[x], 'bjetCSVVeto', 'Number of b-jets', 5, -0.5, 4.5)
#      self.book(names[x], 'muVetoPt5IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
#      self.book(names[x], 'muVetoPt15IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
#      self.book(names[x], 'tauVetoPt20Loose3HitsVtx', 'Number of extra taus', 5, -0.5, 4.5)
#      self.book(names[x], 'eVetoCicTightIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)

#      self.book(names[x], 'jetVeto20', 'Number of extra jets', 5, -0.5, 4.5)
#      self.book(names[x], 'jetVeto30', 'Number of extra jets', 5, -0.5, 4.5)
      
#      self.book(names[x], "mPhiMtPhi", "", 100, 0,4)
#      self.book(names[x], "mPhiMETPhiMVA", "", 100, 0,4)
#      self.book(names[x], "tPhiMETPhiMVA", "", 100, 0,4)
#      self.book(names[x], "mPhiMETPhiType1", "", 100, 0,4)
#      self.book(names[x], "tPhiMETPhiType1", "", 100, 0,4)
    
        
  def correction(self, row):
    return mc_corrector(row)

  def fakeRateMethod(self, row):
    return getFakeRateFactor(row)
      
#  def fill_histo_Ztt(self, row, name='gg', fakeRateN=0, fakeRateD=0):
#    histos = self.histograms
#    weight = self.correction(row)
#    if fakeRateN != 0 and fakeRateD !=0:
#      weight = weight*self.fakeRateMethod(row, fakeRateN, fakeRateD) 
#    histos[name+'/Ztautau_before0'].Fill(row.isZtautau, weight)
                              
  def fill_histos(self, row, name='gg', FRflag=0):
    histos = self.histograms
    weight = self.correction(row)
#    if isEmbed == True:
#      weight = weight*row.EmbPtWeight
    if FRflag==1:
      weight = weight*self.fakeRateMethod(row)
        
    histos[name+'/weight'].Fill(weight)
 #   histos[name+'/weight_nopu'].Fill(self.correction(row))
 #   histos[name+'/rho'].Fill(row.rho, weight)
 #   histos[name+'/nvtx'].Fill(row.nvtx, weight)
 #   histos[name+'/prescale'].Fill(row.doubleMuPrescale, weight)
    
#    histos[name+'/NUP'].Fill(row.NUP, weight)
#    histos[name+'/ZTauTau'].Fill(row.isZtautau, weight)
    
    histos[name+'/mPt'].Fill(row.mPt, weight)
    histos[name+'/mEta'].Fill(row.mEta, weight)
 #   histos[name+'/mMtToMVAMET'].Fill(row.mMtToMVAMET,weight)
 #   histos[name+'/mMtToPfMet_type1'].Fill(row.mMtToPfMet_type1,weight)
 #   histos[name+'/mCharge'].Fill(row.mCharge, weight)
 #   histos[name+'/mRelPFIsoDBHtt'].Fill(row.mRelPFIsoDBHtt, weight)
    
    histos[name+'/tPt'].Fill(row.tPt, weight)
    histos[name+'/tJetPt'].Fill(row.tJetPt, weight)
    histos[name+'/tEta'].Fill(row.tEta, weight)
 #   histos[name+'/tMtToMVAMET'].Fill(row.tMtToMVAMET,weight)
    histos[name+'/tMtToPfMet_type1'].Fill(row.tMtToPfMet_type1,weight)
 #   histos[name+'/tCharge'].Fill(row.tCharge, weight)
    histos[name+'/tDecayMode'].Fill(row.tDecayMode, weight)
#    histos[name+'/tLooseIso'].Fill(row.row.tByLoose###, weight)
    histos[name+'/tByLooseCombinedIsolationDeltaBetaCorr3Hits'].Fill(row.tByLooseCombinedIsolationDeltaBetaCorr3Hits, weight)
    histos[name+'/tByTightCombinedIsolationDeltaBetaCorr3Hits'].Fill(row.tByTightCombinedIsolationDeltaBetaCorr3Hits, weight)
    
#    histos[name+'/LT'].Fill(row.LT,weight)
#    histos[name+'/fullMT_mva'].Fill(fullMT(row.mva_metEt,row.mPt,row.tPt,row.mva_metPhi, row.mPhi, row.tPhi),weight)
#    histos[name+'/fullMT_type1'].Fill(fullMT(row.type1_pfMetEt,row.mPt,row.tPt,row.type1_pfMetPhi, row.mPhi, row.tPhi),weight)
    histos[name+'/collMass_type1'].Fill(collMass_type1(row), weight)
    histos[name+'/collMass'].Fill(collMass(row), weight)

    #    histos[name+'/fullPT_mva'].Fill(fullPT(row.mva_metEt,row.mPt,row.tPt,row.mva_metPhi, row.mPhi, row.tPhi),weight)
#    histos[name+'/fullPT_type1'].Fill(fullPT(row.type1_pfMetEt,row.mPt,row.tPt,row.type1_pfMetPhi, row.mPhi, row.tPhi),weight)
    
    histos[name+'/type1_pfMetEt'].Fill(row.type1_pfMetEt,weight)
    histos[name+'/type1_pfMetPhi'].Fill(row.type1_pfMetPhi,weight)
#    histos[name+'/mva_metEt'].Fill(row.mva_metEt,weight)
    
    histos[name+'/raw_pfMetPhi'].Fill(row.raw_pfMetEt,weight)
    histos[name+'/raw_pfMetPhi'].Fill(row.raw_pfMetPhi,weight)
#    histos[name+'/mva_metPhi'].Fill(row.mva_metPhi,weight)
    
#    histos[name+'/m_t_Mass'].Fill(row.m_t_Mass,weight)
#    histos[name+'/m_t_Pt'].Fill(row.m_t_Pt,weight)
#    histos[name+'/m_t_DR'].Fill(row.m_t_DR,weight)
#    histos[name+'/m_t_DPhi'].Fill(row.m_t_DPhi,weight)
#    histos[name+'/m_t_SS'].Fill(row.m_t_SS,weight)
#    histos[name+'/m_t_ToMETDPhi_type1'].Fill(row.m_t_ToMETDPhi_type1,weight)
    
#    histos[name+'/mPixHits'].Fill(row.mPixHits, weight)
#    histos[name+'/mJetBtag'].Fill(row.mJetBtag, weight)
    
#    histos[name+'/bjetVeto'].Fill(row.bjetVeto, weight)
#    histos[name+'/bjetCSVVeto'].Fill(row.bjetCSVVeto, weight)
#    histos[name+'/muVetoPt5IsoIdVtx'].Fill(row.muVetoPt5IsoIdVtx, weight)
#    histos[name+'/muVetoPt15IsoIdVtx'].Fill(row.muVetoPt15IsoIdVtx, weight)
#    histos[name+'/tauVetoPt20Loose3HitsVtx'].Fill(row.tauVetoPt20Loose3HitsVtx, weight)
#    histos[name+'/eVetoCicTightIso'].Fill(row.eVetoCicTightIso, weight)
#    histos[name+'/jetVeto20'].Fill(row.jetVeto20, weight)
#    histos[name+'/jetVeto30'].Fill(row.jetVeto30, weight)
    
#    histos[name+'/mPhiMtPhi'].Fill(deltaPhi(row.mPhi,row.tPhi),weight)
#    histos[name+'/mPhiMETPhiMVA'].Fill(deltaPhi(row.mPhi,row.mva_metPhi),weight)
#    histos[name+'/tPhiMETPhiMVA'].Fill(deltaPhi(row.tPhi,row.mva_metPhi),weight)
#    histos[name+'/mPhiMETPhiType1'].Fill(deltaPhi(row.mPhi,row.type1_pfMetPhi),weight)
#    histos[name+'/tPhiMETPhiType1'].Fill(deltaPhi(row.tPhi,row.type1_pfMetPhi),weight)


################################################################################
#### SELECTIONS ################################################################
################################################################################

  def presel(self, row):
    if not row.singleIsoMu20Pass:
      return False
    
#    if isZeroJet:
#      if row.NUP != 5:           #CHANGE switched from jetVeto30, ~10% difference
#        return False

#    if isZeroJet:
#      if row.jetVeto30 != 0:
#        return False
    
#    if isEmbed:                 
#      if not row.doubleMuPass:
#        return False
            
          
    return True

  def kinematics(self, row):
    if row.mPt<25 :
      return False
    if row.tPt<20 :
      return False
    if abs(row.tEta) >= 2.3 :
      return False
    if abs(row.mEta) >= 2.1:
      return False
    return True

  def gg0(self,row):
#    if row.jetVeto30!=0:
#      return False
    if row.tPt < 35:
      return False
    if row.mPt < 45 :
      return False
    if row.tMtToPfMet_type1>50:
      return False
#       if row.mMtToPfMet_type1<30:
#            return False
#       if (deltaPhi(row.tPhi,row.pfMetPhi)>0.3):
#          return False
    if (deltaPhi(row.mPhi,row.tPhi)<2.7):    
      return False
#       if (deltaPhi(row.mPhi,row.pfMetPhi)<2.7):
#          return False
    return True

#  def noZtt(self,row):
#    if NoZtautau:
#      if row.isZtautau:
#        return False
#    return True

  def oppositesign(self,row):
    if row.mCharge*row.tCharge!=-1:
      return False
    return True

  def samesign(self,row):
    if row.mCharge!=row.tCharge:
      return False
    return True  
    
  def obj1_id(self, row):
    return bool(row.mPFIDTight) and bool(abs(row.mPVDZ) < 0.2)
  
  def obj2_id(self, row):
    return  row.tAgainstElectronMediumMVA5 and row.tAgainstMuonTight3 and row.tDecayModeFinding

  def vetos(self,row):
    return bool (row.muVetoPt5IsoIdVtx<1) and bool (row.tauVetoPt20Loose3HitsVtx<1)
#and bool (row.eVetoMVAIso) #ALWAYS ZERO
  def obj1_iso(self, row):
    return bool(row.mRelPFIsoDBDefault <0.1)

  def obj2_iso_loose(self, row):
    return  row.tByLooseCombinedIsolationDeltaBetaCorr3Hits
  def obj2_iso_tight(self, row):
    return  row.tByTightCombinedIsolationDeltaBetaCorr3Hits
        
  def obj1_antiiso(self, row):
    return bool(row.mRelPFIsoDBDefault >0.2) 
  
  def obj2_antiiso_tight(self, row):
    return not row.tByTightCombinedIsolationDeltaBetaCorr3Hits


################################################################################
#### RUN!    ###################################################################
################################################################################


  def process(self):
    event =0                 #CHANGE removed this feature to match Aaron
    sel=False
    for row in self.tree:
      if event!=row.evt:   # This is just to ensure we get the (Mu,Tau) with the highest Pt
        event=row.evt    # In principle the code saves all the MU+Tau posibilities, if an event has several combinations
        sel = False      # it will save them all.
      if sel==True:
        continue

      if not self.presel(row):
        continue
      if not self.obj1_id(row):
        continue
      #if not self.obj2_id(row):
        #continue
      if not self.kinematics(row):
        continue
      if not self.vetos(row):
        continue

      sel = True

      #zero jet
      if self.obj1_iso(row) and row.jetVeto30==0:
              
        #Old tau ID and Iso
        if self.obj2_id(row):

          if self.obj2_iso_tight(row) and self.oppositesign(row):
            self.fill_histos(row,'pre_gg')
            if self.gg0(row):
              self.fill_histos(row,'gg')
              
          if self.obj2_iso_loose(row) and not self.obj2_iso_tight(row) and self.oppositesign(row):
            self.fill_histos(row,'pre_gg_looseNotTightTau', 1)
            if self.gg0(row):
              self.fill_histos(row,'gg_looseNotTightTau', 1)
              
          if self.obj2_iso_tight(row) and self.samesign(row):
            self.fill_histos(row,'pre_gg_ss')
            if self.gg0(row):
              self.fill_histos(row,'gg_ss')

          if self.obj2_iso_loose(row) and not self.obj2_iso_tight(row) and self.samesign(row):
            self.fill_histos(row,'pre_gg_ssLooseNotTightTau', 1)                    
            if self.gg0(row):
              self.fill_histos(row,'gg_ssLooseNotTightTau', 1)                    
                
                
  def finish(self):
    self.write_histos()
