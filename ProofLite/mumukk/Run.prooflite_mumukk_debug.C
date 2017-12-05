{

#include "TTree.h"
#include "TDSet.h"
#include "TProof.h"
#include "TString.h"


  // INPUT DATA SAMPLE ON LOCAL DISK

  TDSet* dataset = new TDSet("TTree", "X_data", "mkcands");
  //
  dataset->Add("/lustre/cms/store/user/adiflori/MuOniaParked/runB_resplit_Oct17/runB_split_00.root");

  // dataset->Add("/Users/adrianodiflorio/Documents/Git/X4140/ProofLite/Y4140_testrootuple.root");
  TString selector = "mumukk";
  TProof *p = TProof::Open("workers=1"); // 12 workers for qsub

  // Processing
  cout << ">> Processing " << selector << " ... " << endl;



  TString selectorplus = selector;
  selectorplus += ".C+";
  p->Process(dataset, selectorplus);

}
