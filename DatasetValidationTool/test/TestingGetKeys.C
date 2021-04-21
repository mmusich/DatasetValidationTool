void Recursion(TDirectory* , const char*, const char*); 
float* GetParameters(ROOT::RDataFrame , TString);

void TestingGetKeys(const char* fileName1, const char* fileName2)
//void TestingGetKeys()
{
  std::cout<<"Testing Started ....."<<std::endl;
//  char fileName1[20], fileName2[20];
//  fileName1="Temp.root";  fileName2="Temp.root";
  TString InputFile = Form("%s",fileName1);
  TFile *file=new TFile(InputFile.Data());
  
  TString target1=file->GetPath();
 // check folder exist or not
  gSystem->Exec("mkdir Plots");
  
  file->GetListOfKeys()->Print();
  Recursion(file, fileName1, fileName2); 
}


void Recursion(TDirectory* target, const char* filename1, const char* filename2)
{
  TFile *f1=TFile::Open( filename1, "r" );
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );
  f1->cd(path);
  TDirectory *temp= gDirectory;
  std::cout<<"Checking initial Get Keys Here: "<<temp->GetPath()<<std::endl;
   
  TIter next(temp->GetListOfKeys());
  TKey *key;
  while ((key=(TKey*)next()))
    {
      printf("key: %s points to an object of class: %s \n",
	     key->GetName(),
	     key->GetClassName());

      TObject *obj = key->ReadObj();
      if ( obj->IsA()->InheritsFrom( TDirectory::Class()))
	{
	  std::cout << "Found subdirectory " << obj->GetName() <<std::endl;
	  TDirectory *subdir = (TDirectory*)obj;
	  std::cout<<subdir->GetPath()<<std::endl;
	  Recursion( subdir ,  filename1, filename2);
	}
      else if (obj->IsA()->InheritsFrom( TTree::Class() ))
	{	  
	  auto CSVfileName = "HistParameters.csv";
	  auto df = ROOT::RDF::MakeCsvDataFrame(CSVfileName);
	  
	  std::cout<<"Tree found"<<std::endl;
       	  TTree *tree = (TTree*)obj;
	  TDirectory *TreeDir = tree->GetDirectory();
	  TString DirPath( (char*)strstr( TreeDir->GetPath(), ":" ) );
	  DirPath.Remove( 0, 2 );
	  
          TFile *file2=TFile::Open( filename2, "r" );
	  TString treeName = key->GetName();
	  TString cmd2 = DirPath+"/"+treeName;
          // Check whether tree exist in other files
	  TTree *tree2 = (TTree*)file2->Get(cmd2.Data());		      
	  
	  TH1F *h1,*h2;
	  TCanvas c("c","c",800,600);

	  TIter next2(tree->GetListOfBranches());
	  TKey *key2;
	  while ((key2=(TKey*)next2()))
	    { 
	       
	      TString BranchName=key2->GetName();

	      if(treeName== "Event")
		{
		  float* Par = GetParameters(df, BranchName.Data());
		  int Bins = Par[0];
		  float xMin = Par[1];
		  float xMax = Par[2];	       
		  
		  // loop over all files to draw all branches into a single histogram 
		  TString Drawcmd1 = Form("%s>>hn(%i,%f,%f)",BranchName.Data(),Bins,xMin,xMax);
		  tree->Draw(Drawcmd1.Data(),"","goff",1000,0);
		  h1 = (TH1F*)gDirectory->Get("hn");
		  TString Drawcmd2 = Form("%s>>hm(%i,%f,%f)",BranchName.Data(),Bins,xMin,xMax);
		  tree2->Draw(Drawcmd2.Data(),"","goff",100,0);
		  h2 = (TH1F*)gDirectory->Get("hm");

		  h1->SetLineColor(kBlue);
		  h1->SetLineWidth(3);
		  h1->DrawNormalized();
		  h2->SetLineColor(kRed);
		  h2->SetLineWidth(3);
		  h2->DrawNormalized("Same");
 	     // save in pdf, root and c

		  TString PNGname = treeName + "_" + BranchName + ".png";
		  c.SaveAs(PNGname.Data());
		  TString cmd4 = "mv  " + PNGname + "  Plots";
		  gSystem->Exec(cmd4.Data());
		  
		  c.Clear();    
		  h1->Clear();
		  h2->Clear();	  
		 }	      
	    }
	}
      else
	{std::cout<<"no new directory:"<<std::endl;
	}
      
    }
}

float* GetParameters(ROOT::RDataFrame df, TString var)
{
  static float HistPar[3];  
  TString cmd = "Variable == \""+ var+ "\"";
  std::cout<<"THIS IS CMD DATA CHECK: "<<cmd.Data()<<std::endl;
  auto filteredEvents = df.Filter(cmd.Data());
  auto bins = filteredEvents.Take<Long64_t>("Bins");
  auto xMin = filteredEvents.Take<double>("Xmin");
  auto xMax = filteredEvents.Take<double>("Xmax");

  std::cout<<"Bins value is: "<<bins->at(0)<<std::endl;
  std::cout<<"xMin value is: "<<xMin->at(0)<<std::endl;
  std::cout<<"xMax is: "<<xMax->at(0)<<std::endl;
  HistPar[0]=bins->at(0); 
  HistPar[1]=xMin->at(0);
  HistPar[2]=xMax->at(0);
  bins->clear();
  xMin->clear();
  xMax->clear();
  return HistPar;
  
}
