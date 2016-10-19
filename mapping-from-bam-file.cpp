#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace BamTools;

struct cut_site {
  int start;
  int end;
  int left_len;
  int right_len;
  int lcounts;
  int rcounts;
};

int main(int argv,char **argc) {
  // Read fragments
  std::ifstream fragments(argc[1]);
  std::map<std::string,std::vector<cut_site*> > cuts;
  
  if(!fragments) {
    std::cerr << "Cannot read " << argc[1] << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "Reading fragments" << std::endl;
  
  std::string chr;
  
  while(fragments >> chr) {
    int start, end, lright, lleft;
    std::string lr, ll;
    
    fragments >> start;
    fragments >> end;
    fragments >> lr;
    fragments >> ll;
    
    if(lr == "NA") {
      lright = -1;
    } else {
      std::stringstream ss;
      
      ss << lr;
      ss >> lright;
    }
    
    if(ll == "NA") {
      lleft = -1;
    } else {
      std::stringstream ss;
      
      ss << ll;
      ss >> lleft;
    }
    
    cut_site *s = new cut_site;
    
    s->start = start;
    s->end = end;
    s->left_len = lleft;
    s->right_len = lright;
    s->lcounts = s->rcounts = 0;
    
    if(cuts.count(chr) == 0) {
      cuts[chr] = std::vector<cut_site*>();
    }
    
    cuts[chr].push_back(s);
  }
  
  fragments.close();
  
  std::cout << "Reading BAM file" << std::endl;
  
  BamReader bam;
  
  if(!bam.Open(argc[2])) {
    std::cerr << "Cannot find " << argc[2] << std::endl;
    return EXIT_FAILURE;
  }
  
  BamAlignment align;
  int curchr;
  std::string curchrs;
  int curidx = 0;
  int cutlen = -1;
  int readlen = -1;
  
  while(bam.GetNextAlignment(align)) {
    
    if(!align.IsMapped()) continue;
    
    if(readlen == -1) {
      readlen = align.QueryBases.length();
    }
    
    if(curchr != align.RefID) {
      int refids = bam.GetReferenceCount();
      const RefVector& refs = bam.GetReferenceData();
      
      for( RefData r : refs) {
        int tmpid = bam.GetReferenceID(r.RefName);
        //std::cerr << r.RefName << " " << tmpid << " " << align.RefID << std::endl;
        
        if(tmpid == align.RefID) {
          curchr = tmpid;
          curchrs = r.RefName;
          cutlen = cuts[curchrs].size();
          break;
        }
      }
      
      if(curchr != align.RefID) {
        // sanity check
        std::cerr << "Could not find a reference ID somehow" << std::endl;
        return EXIT_FAILURE;
      }
      
      curidx = 0;
      cutlen = cuts[curchrs].size();
    }
    
    cut_site *s;
    int previdx = curidx;
    bool found = false;
    
    while( curidx < cutlen ) {
      if( cuts[curchrs][curidx]->start == align.Position) {
        cuts[curchrs][curidx]->rcounts++;
        found = true;
        break;
      }
      if(cuts[curchrs][curidx]->end == (align.Position+readlen)) {
        cuts[curchrs][curidx]->lcounts++;
        found = true;
        break;
      }
      
      curidx++;
    }
    
    if(!found) curidx = previdx; // somehow a non-RE read got into the system
  }
  
  std::cout << "Writing output" << std::endl;
  
  std::ofstream raw(argc[3]);
  std::ofstream filtered(argc[4]);
  
  if(!raw ){
    std::cerr << "Cannot write!" << std::endl;
    return EXIT_FAILURE;
  }
  
  if(!filtered) {
    std::cerr << "Cannot write!" << std::endl;
    return EXIT_FAILURE;
  }
  
  raw << "track type=wiggle_0" << std::endl;
  filtered << "track type=wiggle_0" << std::endl;
  
  for(std::map<std::string,std::vector<cut_site*> >::iterator it = cuts.begin(); it != cuts.end(); it++ ) {
    std::string chr = it->first;
    bool foundr = false, foundf = false;
    std::string header = "variableStep chrom=" + chr;
    for( cut_site *s : it->second ) {
      if(s->lcounts > 0 ) {
        if(!foundr) {
          raw << header << std::endl;
          foundr = true;
        }
        raw << s->start << "\t" << s->lcounts << std::endl;
        
        if(s->left_len != -1 ) {
          if(!foundf) {
            filtered << header << std::endl;
            foundf = true;
          }
          filtered << s->start << "\t" << s->lcounts << std::endl;
        }
      }
      
      if(s->rcounts > 0 ) {
        if(!foundr) {
          raw << header << std::endl;
          foundr = true;
        }
        raw << s->end << "\t" << s->rcounts << std::endl;
        
        if(s->right_len != -1 ) {
          if(!foundf) {
            filtered << header << std::endl;
            foundf = true;
          }
          filtered << s->end << "\t" << s->rcounts << std::endl;
        }
      }
    }
  }
  
  raw.close();
  filtered.close();
  
  return 0;
}