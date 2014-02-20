#triggers list
# a text file with a line for each trigger to be displayed in the Muon_trk QA output file with the following form:
#
# name  B-trigger                AC-trigger                E-trigger
# MB7 CINT7-B-NOPF-ALLNOTRD   CINT7-AC-NOPF-ALLNOTRD   CINT7-E-NOPF-ALLNOTRD
# NOSHOW CINT7-I-NOPF-ALLNOTRD   notrigger                notrigger
# MB1 CINT1-B-NOPF-ALLNOTRD   CINT1-AC-NOPF-ALLNOTRD   CINT1-E-NOPF-ALLNOTRD
# MUONLPT CMUS7-B-NOPF-MUON       CMUS7-AC-NOPF-MUON       CMUS7-E-NOPF-MUON
# MUONHPT CMUSH7-B-NOPF-MUON      CMUSH7-AC-NOPF-MUON      CMUSH7-E-NOPF-MUON
# MUONUNLIKE CMUU7-B-NOPF-ALLNOTRD   CMUU7-AC-NOPF-ALLNOTRD   CMUU7-E-NOPF-ALLNOTRD
# MUONLIKE CMUL7-B-NOPF-MUON       CMUL7-AC-NOPF-MUON      CMUL7-E-NOPF-MUON
# MUONLPTS CMUS7-S-NOPF-MUON       CMUS7-ACE-NOPF-MUON       notrigger
# MUONHPTS CMUSH7-S-NOPF-MUON      CMUSH7-ACE-NOPF-MUON      notrigger 
# 

# with passX, the official QA output name is QAresults.root
# with cpassX, the official QA output name is QAresults_barrel.root (for calibration triggers) and QAresults_outer.root (for muon triggers) 

$ALICE_ROOT/PWGPP/MUON/lite/runQA.sh \
-f -o QAresults_outer.root \
-i triggerList.txt \
runList.txt cpass1 \
alien:///alice/data/2012/LHC12d

