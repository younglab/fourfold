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
  
  if(argv<9) {
    std::cerr << "not enough arguments" << std::endl;
    return EXIT_FAILURE;
  }
  
  std::stringstream coords;
  coords << argc[6] << " " << argc[7] << " " << argc[8];
  
  std::string primerchr;
  int primers, primere;
  
  coords >> primerchr;
  coords >> primers;
  coords >> primere;
  
  std::cerr << "test: " << primerchr << " " << primers << " " << primere << std::endl;
  
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
  int totalreads = 0;
  int mappedreads = 0;
  int transreads = 0;
  int cisreads = 0;
  int abnormalreads = 0;
  std::vector<cut_site*> vec;
  
  while(bam.GetNextAlignment(align)) {
    totalreads++;
    
    if(!align.IsMapped()) continue;
    mappedreads++;
    
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
      vec = cuts[curchrs];
      cutlen = vec.size();
    }
    
    cut_site *s;
    int previdx = curidx;
    bool found = false;
    
    while( curidx < cutlen ) {
      if( vec[curidx]->start == align.Position) {
        vec[curidx]->rcounts++;
        found = true;
        break;
      }
      if(vec[curidx]->end == (align.Position+readlen)) {
        vec[curidx]->lcounts++;
        found = true;
        break;
      }
      
      curidx++;
    }

    if(!found) {
      curidx = previdx; // somehow a non-RE read got into the system
      abnormalreads++;
    } else {
      if(curchrs == primerchr) {
        cisreads++;
      } else {
        transreads++;
      }
    }
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
  
  int filteredcis = 0;
  int filteredtrans = 0;
  
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
          if( chr == primerchr ) filteredcis += s->lcounts;
          else filteredtrans += s->lcounts;
          
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
          if( chr == primerchr ) filteredcis += s->rcounts;
          else filteredtrans += s->rcounts;
        }
      }
    }
  }
  
  raw.close();
  filtered.close();
  
  std::ofstream stats(argc[5]);
  
  if(!stats) {
    std::cerr << "Not a valid stats file" << std::endl;
    return EXIT_FAILURE;
  }
  
  stats << "Total Reads: " << totalreads << std::endl;
  stats << "Mapped Reads: " << mappedreads << std::endl;
  stats << "Cis Reads: " << cisreads << std::endl;
  stats << "Trans Reads: " << transreads << std::endl;
  stats << "Abnormal Reads: " << abnormalreads << std::endl;
  stats << "Filtered Cis Reads: " << filteredcis << std::endl;
  stats << "Filtered Trans Reads: "  << filteredtrans << std::endl;
  
  std::vector<cut_site*> sites_in_cis = cuts[primerchr];
  int i;
  
  for( i = 0; i < sites_in_cis.size(); i++ ) {
    cut_site *r = sites_in_cis[i];
    if((r->start+1) == primers || (r->end+1) == primere) break; // right now the cutting coordinates are stored in a 0-offset while the primer coordinates are in a 1-offset, need to synchronize better
  }
  
  if( i == sites_in_cis.size()) {
    stats << "Self-ligation Reads: NA" << std::endl;
    stats << "Non-cut Reads: NA" << std::endl;
  } else {
    int l = i-1;
    int r = i+1;
    
    cut_site *pl = sites_in_cis[l];
    cut_site *p = sites_in_cis[i];
    cut_site *pr = sites_in_cis[r];
    
    if((p->start+1) == primers) {
      stats << "Self-ligation Reads: " << pr->lcounts << std::endl;
      stats << "Non-cut Reads: " << p->lcounts << std::endl;
    } else {
      stats << "Self-ligation Reads: " << pl->rcounts << std::endl;
      stats << "Non-cut Reads: " << p->rcounts << std::endl;
    }
  }
  
  
  stats.close();
  
  return 0;
}