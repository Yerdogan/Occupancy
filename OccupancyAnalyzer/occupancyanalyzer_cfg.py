import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")

#FileService for histograms
process.TFileService = cms.Service("TFileService",
    fileName=cms.string('occ_output.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:0EAE8FA0-0001-E311-AEAD-003048679076.root'
    )
    #, eventsToProcess = cms.untracked.VEventRange('1:63:62021')
)

process.JetGenJetMatch = cms.EDProducer("GenJetMatcher",    # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = cms.InputTag("ak5PFJetsJEC"),       # RECO jets (any View<Jet> is ok)
    matched     = cms.InputTag("ak5GenJets"),        # GEN jets  (must be GenJetCollection)
    mcPdgId     = cms.vint32(),                      # n/a
    mcStatus    = cms.vint32(),                      # n/a
    checkCharge = cms.bool(False),                   # n/a
    maxDeltaR   = cms.double(0.1),                   # Minimum deltaR for the match
    maxDPtRel   = cms.double(3.0),                   # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),
)

process.demo = cms.EDAnalyzer('OccupancyAnalyzer'
)


process.p = cms.Path(process.JetGenJetMatch * process.demo)
