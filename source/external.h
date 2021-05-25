#ifndef __phasex__external__
#define __phasex__external__

#include <iostream>
#include <cstdio>
#include <map>
#include <unordered_map>
#include <vector>

using namespace std;


extern string Program_name;
extern string Program_version;
extern string Program_author;
extern string Program_contrib;
extern string Program_date;
extern string v_website;
extern string v_message;
extern int ready;
extern int debug;

extern int v_quiet;
extern string v_vcf;
extern string homedir;
extern string v_threads;
extern int v_useconfig;
extern int screen_size;
extern string Program_contact;
extern vector <string> warnings;
extern int v_iteration;
extern int v_replicates;
extern float v_threshold;
extern float v_select;
extern string v_beagle;
extern string v_shapeit;
extern string v_output;
extern string v_region;
extern string v_map;;
extern int multiallelic;
extern int onlybiallelic;
extern string v_input;


extern map <pair<string,string>,string> db;
extern unordered_map <string,string> snp;
extern vector <string> sample_list;
extern vector <string> snp_list;
extern map <string,int> used_pairs;
extern map <pair<string,string>,int> final_haplotypes_vcf;
extern vector <string> list_next_round;
extern unordered_map <string,string> fixed_haplotypes; //********
extern int last_pass;
extern vector <string> passed_samples;

extern string shapeit_scheme;
extern string shapeit_useps;
extern string shapeit_pbwt_depth;
extern string shapeit_other_parameters;

extern map <pair<string,string>,int> final_haplotypes;
extern int min_fix;

//extern string v_rbp;
//extern string v_ws;
//extern string v_hc;
//extern string v_beagle;
//extern string v_phase;
//extern string v_beagle_parameters;
//extern string v_phase_parameters;
//extern float v_min_threshold;
//extern int v_cleanup;
//extern int v_debug;
//extern int v_nocombined;

//extern int v_replicates_phase;
//extern string v_output;
//extern float v_select;
//extern string v_memory;

/*
extern map <pair<string,string>,string> db;
extern map <string,string> snp;
extern vector <string> sample_list;
extern vector <string> snp_list;

extern map <pair<string,int>,string> knowndb;
extern map <string,string> knownblocks;
extern map <string,int> randblock;

extern map <pair<string,string>,int> final_haplotypes_rbp;
extern map <pair<string,string>,int> final_haplotypes_ws;
extern map <pair<string,string>,int> final_haplotypes_vcf;
extern map <pair<string,string>,int> final_haplotypes_combined;
extern map <string,string> haplotype_block_hc;
extern map <string,int> used_pairs;
extern int last_pass;
extern vector <string> list_next_round;
extern map <string,string> fixed_haplotypes;
extern int v_unphase;

extern string v_input;
*/
#endif /* defined(__phasex__external__) */

