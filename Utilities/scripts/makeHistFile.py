
#!/usr/bin/env python
import ROOT

## from the directory python, include these files
from python import SelectorTools
from python import UserInput
from python import OutputTools
from python import ConfigureJobs
from python import HistTools

#additional python libraries
import os
import sys
 


## get command line arguments: sets parameters based off of what is
## given to initiate run of this program
def getComLineArgs():
    parser = UserInput.getDefaultParser()
    parser.add_argument("--proof", "-p", 
        action='store_true', help="Don't use proof")
    parser.add_argument("--lumi", "-l", type=float,
        default=35.87, help="luminosity value (in fb-1)")
    parser.add_argument("--output_file", "-o", type=str,
        default="test.root", help="Output file name")
    parser.add_argument("--test", action='store_true',
        help="Run test job (no background estimate)")
    parser.add_argument("--noHistConfig", action='store_true',
        help="Don't rely on config file to specify hist info")
    parser.add_argument("--output_selection", type=str,
        default="", help="Selection stage of output file "
        "(Same as input if not give)")
    parser.add_argument("-b", "--hist_names", 
                        type=lambda x : [i.strip() for i in x.split(',')],
                        default=["all"], help="List of histograms, "
                        "as defined in AnalysisDatasetManager, separated "
                        "by commas")
    return vars(parser.parse_args())

## anything you do with ROOT will not require permission to 
## access the display and you can generate plots, save, etc. 
ROOT.gROOT.SetBatch(True)

# calling the function getComLineArgs() given above 
args = getComLineArgs()

#calling on the function getManagerPath() that is in the ConfigureJobs file
manager_path = ConfigureJobs.getManagerPath()

#giving a temporary file name of output_file
tmpFileName = args['output_file']

#opens or creates a local ROOT file calling it fOut. The recreate option creates a new file; 
# if the file already exists, it will be overwritten
fOut = ROOT.TFile(tmpFileName, "recreate")


#same as above; but no specified mode in which the file is opened
fScales = ROOT.TFile('data/scaleFactors.root')
mCBTightFakeRate = fScales.Get("mCBTightFakeRate")
eCBTightFakeRate = fScales.Get("eCBTightFakeRate")

# not a clue......  #################################################################
useSvenjasFRs = False
useJakobsFRs = False
if useSvenjasFRs:
    mCBTightFakeRate = fScales.Get("mCBTightFakeRate_Svenja")
    eCBTightFakeRate = fScales.Get("eCBTightFakeRate_Svenja")
elif useJakobsFRs:
    mCBTightFakeRate = fScales.Get("mCBTightFakeRate_Jakob")
    eCBTightFakeRate = fScales.Get("eCBTightFakeRate_Jakob")

# For medium muons
#// mCBMedFakeRate.SetName("fakeRate_allMu") // 
if mCBTightFakeRate:
    mCBTightFakeRate.SetName("fakeRate_allMu")
if eCBTightFakeRate:
    eCBTightFakeRate.SetName("fakeRate_allE")

# getting the value for the specified key (in quotes) from scaleFactor file; 
# reminder: this is how it works->  dictionary = {'key': 'value', etc.}
muonIsoSF = fScales.Get('muonIsoSF')
muonIdSF = fScales.Get('muonTightIdSF')
electronTightIdSF = fScales.Get('electronTightIdSF')
electronGsfSF = fScales.Get('electronGsfSF')
pileupSF = fScales.Get('pileupSF')

# obtaining another root file calling it fPrefireEfficiency 
# and next getting the value for the specified key (like above)
#//fPrefireEfficiency = ROOT.TFile('data/Map_Jet_L1FinOReff_bxm1_looseJet_JetHT_Run2016B-H.root')//
fPrefireEfficiency = ROOT.TFile('data/Map_Jet_L1FinOReff_bxm1_looseJet_SingleMuon_Run2016B-H.root')
prefireEff = fPrefireEfficiency.Get('prefireEfficiencyMap')

# creating a list of values to be gotten later on
fr_inputs = [eCBTightFakeRate, mCBTightFakeRate,]
sf_inputs = [electronTightIdSF, electronGsfSF, muonIsoSF, muonIdSF, pileupSF, prefireEff]

# we are looking at what we entered into the command line (see what args is near the top of this file)
# if we entered output_selection, call it now args['selection']; returns the value in the first index
# which is split based on _. i.e. a_b will now be a b
if args['output_selection'] == '':
    args['output_selection'] = args['selection']
selection = args['output_selection'].split("_")[0]

# if that value is Inclusive2Jet.....
if selection == "Inclusive2Jet":
    selection = "Wselection"
    print "Info: Using Wselection for hist defintions"

# '/'.join joins the values in args['analysis'], selection with a slash between them...assuming some sort of file path? 
# in order to evaluate next line which calls the User Input file, getHistInfo function
analysis = "/".join([args['analysis'], selection])
hists, hist_inputs = UserInput.getHistInfo(analysis, args['hist_names'], args['noHistConfig'])

# ROOT.TNamed contains the essential elements to identify a derived object; in this case 
# the name is "selection" and the title is args['output_selection']
tselection = [ROOT.TNamed("selection", args['output_selection'])]

#creating the channels array
nanoAOD = True
channels = ["Inclusive"] if nanoAOD else ["eee", "eem", "emm", "mmm"]

# again seeing if command line has -p; if so, it is given as TRUE. Then, 
# we create a proof environment -> stars a set of worker servers 
# not sure why.... ################################################################
if args['proof']:
    ROOT.TProof.Open('workers=12')

# if we have our data, with no? FakeRate and not completing a test run? 
# we look at the background. Calling functions in the file Configure Jobs
if "WZxsec2016" in analysis and "FakeRate" not in args['output_selection'] and not args['test']:
    background = SelectorTools.applySelector(["WZxsec2016data"] +
        ConfigureJobs.getListOfEWKFilenames() + ["wz3lnu-powheg"] +
        ConfigureJobs.getListOfNonpromptFilenames(), 
            "WZBackgroundSelector", args['selection'], fOut, 
            extra_inputs=sf_inputs+fr_inputs+hist_inputs+tselection, 
            channels=channels,
            addSumweights=False,
            nanoAOD=nanoAOD,
            proof=args['proof'])

# setting up a dictionary labeled selector map
selector_map = {
    "WZxsec2016" : "WZSelector",
    "Zstudy" : "ZSelector",
}

#using function applySelector defined in Selector Tools file. 
mc = SelectorTools.applySelector(args['filenames'], selector_map[args['analysis']], 
        args['selection'], fOut, 
        extra_inputs=sf_inputs+hist_inputs+tselection, 
        channels=channels,
        nanoAOD=nanoAOD,
        addSumweights=True, proof=args['proof']) 

# if we said test, we are good to go and we close the output file and exit the code
if args['test']:
    fOut.Close()
    sys.exit(0)

# using function makeCompositeHists defined in HistTools file; we then write this onto the file ____ 
# and promptly close the file?? alldata
alldata = HistTools.makeCompositeHists(fOut,"AllData", 
    ConfigureJobs.getListOfFilesWithXSec(["WZxsec2016data"], manager_path), args['lumi'],
    underflow=False, overflow=False)
OutputTools.writeOutputListItem(alldata, fOut)
alldata.Delete()

# similar procedure as defined above for nonpromptmc, ewkmc, and ewkcorr. 
nonpromptmc = HistTools.makeCompositeHists(fOut, "NonpromptMC", ConfigureJobs.getListOfFilesWithXSec( 
    ConfigureJobs.getListOfNonpromptFilenames(), manager_path), args['lumi'],
    underflow=False, overflow=False)
nonpromptmc.Delete()

OutputTools.writeOutputListItem(nonpromptmc, fOut)
ewkmc = HistTools.makeCompositeHists(fOut,"AllEWK", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfEWKFilenames(), manager_path), args['lumi'],
    underflow=False, overflow=False)
OutputTools.writeOutputListItem(ewkmc, fOut)
ewkmc.Delete()

ewkcorr = HistTools.getDifference(fOut, "DataEWKCorrected", "AllData", "AllEWK")
OutputTools.writeOutputListItem(ewkcorr, fOut)
ewkcorr.Delete()
