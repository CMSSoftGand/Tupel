    #include "simpleReader.h"
    void simpleReader(TString fin,TString fout){
      TFile *f = TFile::Open(fin);
      if (f->IsZombie()) {
      printf("Input root files doesn't open, please have a look:\n");
      return;}
    cout<<"this is fileIn:  "<<fin<<"   ;  "<<fout<<endl;
    TTree *t = (TTree*)f->Get("sync");
//    TFile *theFile = new TFile (fout+".root","RECREATE");
//    theFile->cd();  
    ofstream outputFile_Objects;
    outputFile_Objects.open("Object.txt",std::ios::trunc);
      branchAdd(t);
      Int_t nentries(t->GetEntriesFast());

        for (int jentry=0; jentry < nentries; jentry++)
        {
        t->GetEntry(jentry);
        if(jentry%1000==0)cout<<" << "<<jentry<<"/"<<nentries<<endl;

/*      UInt_t Run            = Run;
      ULong64_t LumiSection = LumiSection;
      ULong64_t Event       = Event;
      Bool_t   RecoSuccess  = RecoSuccess;
      Float_t       MassTT  = MassTT;
*/    
if (RecoSuccess)
    //outputFile_Objects<<Run<<":"<< LumiSection << ":"<< Event <<"     MassTT:   "<<MassTT<<"\n"<<endl; 
    outputFile_Objects<<Run<<":"<< LumiSection << ":"<< Event <<endl; 
              
              
              
        }//entries loop
//        theFile->Write();
//        theFile->Close();
        }//function loop
