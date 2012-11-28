import os

################################################################################
#### Fitted fake rate functions ################################################
################################################################################

# Get fitted fake rate functions
frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')
highpt_mu_fr = build_roofunctor(
    #frfit_dir + '/m_wjets_pt20_pfidiso02_muonJetPt.root',
    frfit_dir + '/m_wjets_pt20_h2taucuts_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)
lowpt_mu_fr = build_roofunctor(
    #frfit_dir + '/m_wjets_pt10_pfidiso02_muonJetPt.root',
    frfit_dir + '/m_wjets_pt10_h2taucuts_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)
tau_fr = build_roofunctor(
    frfit_dir + '/t_ztt_pt20_mvaloose_tauPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

highpt_mu_qcd_fr = build_roofunctor(
    frfit_dir + '/m_qcd_pt20_h2taucuts_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)
lowpt_mu_qcd_fr = build_roofunctor(
    frfit_dir + '/m_qcd_pt10_h2taucuts_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

lowpt_e_qcd_fr = build_roofunctor(
    #frfit_dir + '/e_qcd_pt10_mvaidiso03_eJetPt.root',
    frfit_dir + '/e_qcd_pt10_h2taucuts_eJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

tau_qcd_fr = tau_fr ## build_roofunctor(
##     frfit_dir + '/t_ztt_pt20_mvaloose_tauPt.root',
##     'fit_efficiency', # workspace name
##     'efficiency'
## )

# Get 2D fake rates

fr_data_views = data_views.data_views(
    glob.glob(os.path.join('results', os.environ['jobid'], 'FakeRatesMM', '*.root')),
    glob.glob(os.path.join('inputs', os.environ['jobid'], '*.sum')),
)

def get_view(sample_pattern):
    for sample, sample_info in fr_data_views.iteritems():
        if fnmatch.fnmatch(sample, sample_pattern):
            return sample_info['view']
    raise KeyError("I can't find a view that matches %s, I have: %s" % (
        sample_pattern, " ".join(fr_data_views.keys())))

# FR data, subtracting WZ and ZZ.
mu_fr_ewk_2d = TwoDimFakeRate(
    'wjets/pt10/h2taucuts/muonJetVsLeptonPt', 'wjets/pt10/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

mu_fr_qcd_2d = TwoDimFakeRate(
    'qcd/pt10/h2taucuts/muonJetVsLeptonPt', 'qcd/pt10/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

# eta dependent jet-pt vs pt
mu_fr_ewk_2d_f = TwoDimFakeRate(
    'wjets/pt10f/h2taucuts/muonJetVsLeptonPt', 'wjets/pt10f/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))
mu_fr_qcd_2d_f = TwoDimFakeRate(
    'qcd/pt10f/h2taucuts/muonJetVsLeptonPt', 'qcd/pt10f/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

mu_fr_ewk_2d_t = TwoDimFakeRate(
    'wjets/pt10t/h2taucuts/muonJetVsLeptonPt', 'wjets/pt10t/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))
mu_fr_qcd_2d_t = TwoDimFakeRate(
    'qcd/pt10t/h2taucuts/muonJetVsLeptonPt', 'qcd/pt10t/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

mu_fr_ewk_2d_b = TwoDimFakeRate(
    'wjets/pt10b/h2taucuts/muonJetVsLeptonPt', 'wjets/pt10b/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))
mu_fr_qcd_2d_b = TwoDimFakeRate(
    'qcd/pt10b/h2taucuts/muonJetVsLeptonPt', 'qcd/pt10b/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

lowpt_e_fr = build_roofunctor(
    #frfit_dir + '/e_wjets_pt10_mvaidiso03_eJetPt.root',
    frfit_dir + '/e_wjets_pt10_h2taucuts_eJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

