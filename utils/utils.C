#ifndef UTILS
#define UTILS

bool DirectoryExists( const char* pzPath )
{//Check if a directory exists
  if ( pzPath == NULL) return false;
    DIR *pDir;
    bool bExists = false;
    pDir = opendir (pzPath);
    if (pDir != NULL)
    {
        bExists = true;    
        (void) closedir (pDir);
    }
    return bExists;
}


vector <int> get_runs_from_file(string filename)                                                                                                                                                         
{                                                                                                                                                                                                        
  // Open run files                                                                                                                                                                                      
  std::ifstream file(filename.c_str());                                                                                                                                                                  
  // Check if the file opened successfully                                                                                                                                                               
  if (!file.is_open()) {                                                                                                                                                                                 
    std::cerr << "Error opening the file." << std::endl;                                                                                                                                                 
    return {0};                                                                                                                                                                                          
  }                                                                                                                                                                                                      
  string line;                                                                                                                                                                                           
  getline(file, line); // Read the first line                                                                                                                                                            
  file.close();                                                                                                                                                                                          
  vector<int> numbers;                                                                                                                                                                                   
  istringstream iss(line);                                                                                                                                                                               
  string token;                                                                                                                                                                                         
  // Parse the line by comma and store numbers in vector                                                                                                                                                 
  while (getline(iss, token, ',')) {                                                                                                                                                                     
    // Convert string to integer                                                                                                                                                                         
    int num;                                                                                                                                                                                             
    istringstream(token) >> num;                                                                                                                                                                         
    numbers.push_back(num);                                                                                                                                                                              
  }                                                                                                                                                                                                      
  return numbers;                                                                                                                                                                                        
}

#endif
