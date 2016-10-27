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
#include <cmath>

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
  
  if(argv<12) {
    std::cerr << "not enough arguments" << std::endl;
    return EXIT_FAILURE;
  }
  
  const char *fragmentfile = argc[1];
  const char *bamfile = argc[2];
  const char *rawwig = argc[3];
  const char *filteredwig = argc[4];
  const char *countrawtable = argc[5];
  const char *countfilteredtable = argc[6];
  const char *rpmrawtable = argc[7];
  const char *rpmfilteredtable = argc[8];
  const char *statsfile = argc[9];
  const char *pchr = argc[10];
  const char *pstart = argc[11];
  const char *pend = argc[12];
  
  std::stringstream coords;
  coords << pchr << " " << pstart << " " << pend;
  
  std::string primerchr;
  int primers, primere;
  
  coords >> primerchr;
  coords >> primers;
  coords >> primere;
  
  // Read fragments
  std::ifstream fragments(fragmentfile);
  std::map<std::string,std::vector<cut_site*> > cuts;
  
  if(!fragments) {
    std::cerr << "Cannot read " << fragmentfile << std::endl;
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
  
  if(!bam.Open(bamfile)) {
    std::cerr << "Cannot find " << bamfile << std::endl;
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
  
  std::ofstream raw(rawwig);
  std::ofstream filtered(filteredwig);
  std::ofstream crtable(countrawtable);
  std::ofstream cftable(countfilteredtable);
  std::ofstream rrtable(rpmrawtable);
  std::ofstream rftable(rpmfilteredtable);
  
  if(!raw ){
    std::cerr << "Cannot write to " << rawwig << std::endl;
    return EXIT_FAILURE;
  }
  
  if(!filtered) {
    std::cerr << "Cannot write to " << filteredwig << std::endl;
    return EXIT_FAILURE;
  }
  
  if(!crtable) {
    std::cerr << "Cannot write to " << countrawtable << std::endl;
    return EXIT_FAILURE;
  }
  
  if(!cftable) {
    std::cerr << "Cannot write to " << countfilteredtable << std::endl;
    return EXIT_FAILURE;
  }
  
  
  if(!rrtable) {
    std::cerr << "Cannot write to " << rpmrawtable << std::endl;
    return EXIT_FAILURE;
  }
  
  if(!rftable) {
    std::cerr << "Cannot write to " << rpmfilteredtable << std::endl;
    return EXIT_FAILURE;
  }
  
  int filteredcis = 0;
  int filteredtrans = 0;
  
  int coverage200kb = 0;
  int reads200kb = 0;
  int fragends = 0;
  int filteredcoverage200kb = 0;
  int filteredreads200kb = 0;
  int filteredfragends = 0;
  
  raw << "track type=wiggle_0" << std::endl;
  filtered << "track type=wiggle_0" << std::endl;
  
  for(std::map<std::string,std::vector<cut_site*> >::iterator it = cuts.begin(); it != cuts.end(); it++ ) {
    std::string chr = it->first;
    bool foundr = false, foundf = false;
    std::string header = "variableStep chrom=" + chr;
    for( cut_site *s : it->second ) {
      bool close = false;
      
      if(chr == primerchr && (std::abs(s->start-primers) < 100000 || std::abs(s->end-primers) < 100000)) {
        close = true;
        fragends+=2;
        if(s->left_len != -1 ) filteredfragends++;
        if(s->right_len != -1) filteredfragends++;
      }
      
      if(s->lcounts > 0 ) {
        if(!foundr) {
          raw << header << std::endl;
          foundr = true;
        }
        raw << s->start << "\t" << s->lcounts << std::endl;
        
        if(close) {
          coverage200kb++;
          reads200kb+=s->lcounts;
        }
        
        if(s->left_len != -1 ) {
          if(!foundf) {
            filtered << header << std::endl;
            foundf = true;
          }
          filtered << s->start << "\t" << s->lcounts << std::endl;
          if(close) {
            filteredcoverage200kb++;
            filteredreads200kb+=s->lcounts;
          }
          if( chr == primerchr ) filteredcis += s->lcounts;
          else filteredtrans += s->lcounts;
          
          cftable << chr << "\t" << s->start << "\t" << s->lcounts << std::endl;
          rftable << chr << "\t" << s->start << "\t" << double(s->lcounts)/mappedreads*1e6 << std::endl;
          
          
        }
        
        crtable << chr << "\t" << s->start << "\t" << s->lcounts << std::endl;
        rrtable << chr << "\t" << s->start << "\t" << double(s->lcounts)/mappedreads*1e6 << std::endl;
        
      }
      
      if(s->rcounts > 0 ) {
        if(!foundr) {
          raw << header << std::endl;
          foundr = true;
        }
        raw << s->end << "\t" << s->rcounts << std::endl;
        
        
        if(close) {
          coverage200kb++;
          reads200kb+=s->rcounts;
        }
        
        if(s->right_len != -1 ) {
          if(!foundf) {
            filtered << header << std::endl;
            foundf = true;
          }
          filtered << s->end << "\t" << s->rcounts << std::endl;
          
          if(close) {
            filteredcoverage200kb++;
            filteredreads200kb+=s->rcounts;
          }
          if( chr == primerchr ) filteredcis += s->rcounts;
          else filteredtrans += s->rcounts;
          
          cftable << chr << "\t" << s->end << "\t" << s->rcounts << std::endl;
          rftable << chr << "\t" << s->end << "\t" << double(s->rcounts)/mappedreads*1e6 << std::endl;
          
          
        }
        
        crtable << chr << "\t" << s->end << "\t" << s->rcounts << std::endl;
        rrtable << chr << "\t" << s->end << "\t" << double(s->rcounts)/mappedreads*1e6 << std::endl;
        
      }
    }
  }
  
  raw.close();
  filtered.close();
  crtable.close();
  cftable.close();
  rrtable.close();
  rftable.close();
  
  std::ofstream stats(statsfile,std::ofstream::out|std::ofstream::app);
  
  if(!stats) {
    std::cerr << "Cannot write to " << statsfile << std::endl;
    return EXIT_FAILURE;
  }
  
  stats << "Total Reads: " << totalreads << std::endl;
  stats << "Mapped Reads: " << mappedreads << std::endl;
  stats << "Cis Reads: " << cisreads << std::endl;
  stats << "Trans Reads: " << transreads << std::endl;
  stats << "Abnormal Reads: " << abnormalreads << std::endl;
  stats << "Filtered Cis Reads: " << filteredcis << std::endl;
  stats << "Filtered Trans Reads: "  << filteredtrans << std::endl;
  stats << "Reads within 200kb of read primer: " << reads200kb << std::endl;
  stats << "Filtered reads within 200kb of read primer: " << filteredreads200kb << std::endl;
  stats << "Coverage % within 200kb of read primer: " << double(coverage200kb)/fragends*100 << std::endl;
  stats << "Filtered coverage % within 200kb of read primer: " << double(filteredcoverage200kb)/filteredfragends*100 << std::endl;
  
  
  
  std::vector<cut_site*> sites_in_cis = cuts[primerchr];
  int i;
  
  for( i = 0; i < sites_in_cis.size(); i++ ) {
    cut_site *r = sites_in_cis[i];
    if((r->start+1) == primers || (r->end) == primere) break; // right now the cutting coordinates are stored in a 0-offset while the primer coordinates are in a 1-offset, need to synchronize better
  }
  
  int selfligation = -1;
  int noncut = -1;
  
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
      selfligation = pr->lcounts;
      noncut = p->lcounts;
      
      stats << "Self-ligation Reads: " << selfligation << std::endl;
      stats << "Non-cut Reads: " << noncut << std::endl;
    } else {
      selfligation = pl->rcounts;
      noncut = p->rcounts;
      
      stats << "Self-ligation Reads: " << selfligation << std::endl;
      stats << "Non-cut Reads: " << noncut << std::endl;
    }
  }
  
  stats << "Cis % (raw): " << double(cisreads)/(cisreads+transreads)*100.0 << std::endl;
  if( selfligation > -1 && noncut > -1 ) {
    stats << "Cis % (raw, no self or noncut): " << double(cisreads-selfligation-noncut)/(cisreads+transreads-selfligation-noncut)*100 << std::endl;
  } else {
    stats << "Cis % (raw, no self or noncut): NA" << std::endl;
  }
  stats << "Cis % (filtered): " << double(filteredcis)/(filteredcis+filteredtrans)*100.0 << std::endl;
  if( selfligation > -1 && noncut > -1 ) {
    stats << "Cis % (filtered, no self or noncut): " << double(cisreads-selfligation-noncut)/(cisreads+transreads-selfligation-noncut)*100 << std::endl;
  } else {
    stats << "Cis % (filtered, no self or noncut): NA" << std::endl;
  }
  
  
  
  stats.close();
  
  return 0;
}