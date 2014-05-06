void checkTrees()
{
  TFile* f_ = TFile::Open("data/Z_2e/Ntuple_1.root");
  TList* keys_ = f_->GetListOfKeys();
  for( int ii=0; ii<keys_->GetEntries(); ii++ )
    {
      TKey* key_ = (TKey*)keys_->At(ii);
      TObject* obj_ = key_->ReadObj();      
      if( obj_==0 ) continue;
      if(  obj_->IsA()->GetName()!=TString("TTree") ) continue;
      //cout <<  obj_->IsA()->GetName() << endl;
      TTree* tree_ = (TTree*) obj_;
      size_t n_ =  tree_->GetEntries();
      //      if( n_==0 ) continue; 
      TString name_( tree_->GetName() );
      if( !name_.Contains("MET") ) continue;
      cout <<  name_ << " n=" << n_ <<  endl;
      TObjArray* arr_ = tree_->GetListOfBranches();
      int nb_ = arr_->GetEntries();
      for( int jj=0; jj<nb_; jj++ )
	{
	  TBranch* br_ = (TBranch*)arr_->At(jj);
	  cout << "--> " << br_->GetTitle() << endl;
	}
    }
}
