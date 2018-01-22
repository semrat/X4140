{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "xTree", "rootuple");
  //
  dataset->Add("/lustre/cms/store/user/adiflori/MuOnia/phiJpsiTriggersBCDEF.root");

  // dataset->Add("/Users/adrianodiflorio/Documents/Git/X4140/ProofLite/Y4140_testrootuple.root");
  TString selector = "mumumumu";
  TProof *p = TProof::Open("workers=1"); // 12 workers for qsub

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;



  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
