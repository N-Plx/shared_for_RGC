#include "hipo4/reader.h"

void FCup_reading(string path_to_gmn, int run_min, int run_max)
{  
  //Input file
  TChain *chain = new TChain("", "");
  for(int i=run_min; i<=run_max; i++)
    { 
      string filename = path_to_gmn + "gmn_0"+ to_string(i) + ".hipo";
      ifstream my_file(filename);
      if (my_file) {
	cout << "Opening file for run " << i << endl;
	chain->AddFile(filename.c_str());
      } 
      else 
	{
	  cout << "file " + filename + " does not exist" << endl;
	  continue;
	}
    }
    
  //Reading the input
  TObjArray *files = chain->GetListOfFiles();
  cout << "size of files   " << files->GetLast() + 1 << endl;
  
  //Define hipo variables
  hipo::reader reader;
  //reader.setTags(1);
  hipo::dictionary factory;
  hipo::structure particles;
  hipo::event event;
  
  //Useful variables to be read from banks
  int RunNumber, Helicity;
  double FCup_hel_p, FCup_hel_n, FCup_hel_0, FCup_run;
  //just a counter for debugging and all  
  int reading=0;

  //Output 
  ofstream outfile;
  outfile.open ("FCup_"+to_string(run_min)+"_"+to_string(run_max)+".txt");
  outfile << "RunNumber,FCup_pos,FCup_neg,FCup_zero,FCup_run" << endl;

  //Start loop on files
  for (int noffiles = 0; noffiles < files->GetLast() + 1; noffiles++)
    {
      double FCup_run_min=9999, FCup_run_max=-9999;
      //Handling hipo files
      reader.open(files->At(noffiles)->GetTitle());

      //Run scaler beam charge
      clas12reader c12(files->At(noffiles)->GetTitle(),{0});
      c12.scalerReader();//must call this first                                                                                                                                                          
      
      reader.readDictionary(factory);

      hipo::bank HEL_SCALER(factory.getSchema("HEL::scaler"));
      hipo::bank RUN_SCALER(factory.getSchema("RUN::scaler"));
      hipo::bank CONF(factory.getSchema("RUN::config"));
      hipo::bank HEL(factory.getSchema("REC::Event"));

      cout << "File number " << noffiles << " opened" << endl;
      
      //FCup counts from scalers
      FCup_hel_n = 0;
      FCup_hel_p = 0;      
      FCup_hel_0 = 0;
      FCup_run = 0;

      //Start loop on events
      while (reader.next() == true)
        {
      	  //Get hipo info
	  reader.read(event);  
	  event.getStructure(CONF);
	  event.getStructure(HEL_SCALER);
	  event.getStructure(RUN_SCALER);
	  event.getStructure(HEL);

	  reading++;

	  //Event info                                                                                                                                                                                   
          RunNumber = CONF.getInt("run", 0);
	  Helicity = HEL.getInt("helicity", 0);
	  
	  //Scalers information
          for(int row=0; row<HEL_SCALER.getRows(); row++)
	    {
	      int hel = HEL_SCALER.getByte("helicity",row);
	      if(hel == -1) FCup_hel_n+=HEL_SCALER.getFloat("fcupgated",row);
	      if(hel == 1) FCup_hel_p+=HEL_SCALER.getFloat("fcupgated",row);
	      if(hel == 0) FCup_hel_0+=HEL_SCALER.getFloat("fcupgated",row);
	    }

	  if(RUN_SCALER.getRows()>0)
            {
	      float entry = RUN_SCALER.getFloat("fcupgated",0);
	      if(entry>FCup_run_max) FCup_run_max=entry;
	      if(entry<FCup_run_min) FCup_run_min=entry;
	      FCup_run =  FCup_run_max - FCup_run_min;
	    }
	}//END LOOP ON EVENTS
      cout << RUN_SCALER.getFloat("fcupgated",0) << " " << FCup_run_max << " " << FCup_run_min << " " << FCup_run << " " << FCup_hel_n << " " << FCup_hel_p << " " << FCup_hel_0 << endl;
      outfile << RunNumber << "," << FCup_hel_p << "," << FCup_hel_n << "," << FCup_hel_0 << "," << FCup_run << endl;
      cout << "Number of events read : " << reading << endl;
    }//END LOOP ON FILES
  outfile.close();
}



  
