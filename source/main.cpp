#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <map>
#include <unordered_map>
#include <sys/stat.h>
#include <thread>
#include <pwd.h>
#include <algorithm>


#include "functions.hpp"
#include "vcf.hpp"
#include "external.h"
#include "phase-ps.hpp"
#include "recode.hpp"

using namespace std;

auto clock_start = std::chrono::steady_clock::now();
int v_concurentThreadsSupported = std::thread::hardware_concurrency() / 2;
string v_threads = std::to_string(v_concurentThreadsSupported);


string Program_name = "PHASEX";
//string Program_version = "0.7.5"; //support to shapeit 4.1.2 and absense of map=
//string Program_version = "0.7.6"; //support to min_fix
//string Program_version = "0.8"; //support to splited multiallelic
//string Program_version = "0.8.1"; //write vcf using threads
string Program_version = "0.8.2"; //correct splitted multiallelic with different reference alleles
string Program_author = "Erick C. Castelli";
string Program_contact = "erick.castelli@unesp.br";
string Program_contrib = "";
string Program_date = "Nov 27th 2020";
string v_website = "www.castelli-lab.net/apps/phasex";
string v_message = "";
string configfile = "";
int ready = 0;


int v_quiet = 0;
string v_vcf = "";
string homedir = "";
int v_useconfig = 0;
int screen_size = 80;
vector <string> warnings;
int v_iteration = 50;
int v_replicates = 50;
float v_threshold = 0.95;
float v_select = 0.70;
string v_beagle = "";
string v_shapeit = "";
string v_output = "";
string v_region = "";
string v_map = "";
int multiallelic = 0;
int onlybiallelic = 0;
string v_input = "";
int debug = 0;
int min_fix = 2;


string shapeit_scheme = "15b,1p,1b,1p,1b,1p,1b,1p,1b,1p,1b,1p,15m";
string shapeit_useps = "0.0001";
string shapeit_pbwt_depth = "8";
string shapeit_other_parameters = "";

map <pair<string,string>,string> db;
unordered_map <string,string> snp;
vector <string> sample_list;
vector <string> snp_list;
map <string,int> used_pairs;
map <pair<string,string>,int> final_haplotypes_vcf;
vector <string> list_next_round;
unordered_map <string,string> fixed_haplotypes;  //*************************
int last_pass = 0;
vector <string> passed_samples;

map <pair<string,string>,int> final_haplotypes;

//int v_cleanup = 0;
//int v_debug = 0;
//int v_nocombined = 0;
//int v_replicates_phase = 10;
//string v_rbp = "";
//string v_ws = "";
//string v_hc = "";

//string v_phase = "";
//string v_beagle_parameters = "niterations=10";
//string v_phase_parameters = "100 1 100";


//float v_min_threshold = 0.80;

//string v_memory = "12";

//int v_unphase = 0;

/*


map <pair<string,int>,string> knowndb;
map <string,string> knownblocks;
map <string,int> randblock;
map <pair<string,string>,int> final_haplotypes_rbp;
map <pair<string,string>,int> final_haplotypes_ws;

map <pair<string,string>,int> final_haplotypes_combined;




map <string,string> haplotype_block_hc;


*/

void load_config (void)
{
	if (v_useconfig == 0) {
		ifstream config(configfile.c_str());
		if (config)
		{
			for( std::string line; getline( config, line ); )
			{
				
				line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
				vector<string> v_temp;
				boost::split(v_temp,line,boost::is_any_of("="));
				
				if (v_temp[0] == "beagle") {v_beagle = v_temp[1];}
				if (v_temp[0] == "shapeit") {v_shapeit = v_temp[1];}
				if (v_temp[0] == "map") {v_map = v_temp[1];}
			}
		}
		config.close();
	}
	return;
}








void main_help(void)
{
	screen_message (screen_size, 0, "", 1, 0);
	screen_message (screen_size, 0, "", 1, 0);
	screen_message (screen_size, 0, Program_name, 1, 0);
	screen_message (screen_size, 0, "", 1, 0);
	screen_message (screen_size, 2, "* Author  : " + Program_author, 1, 0);
	screen_message (screen_size, 2, "* Contact : " + Program_contact, 1, 0);
	screen_message (screen_size, 2, "* Version : " + Program_version, 1, 0);
	screen_message (screen_size, 2, "* Website : " + v_website, 1, 0);
	screen_message (screen_size, 0, "", 1, 0);
	screen_message (screen_size, 2, "Tools", 1, 0);
	screen_message (screen_size, 5, "phase-ps   Phase variants considering Phase Sets (PS)", 1, 0);
	screen_message (screen_size, 5, "recreate   Recreate the final VCF", 1, 0);
	screen_message (screen_size, 5, "hp-ps      Recode GATK ReadBackedPhasing to PS format", 1, 0);
	screen_message (screen_size, 0, "", 1, 0);
}









int main(int argc, const char * argv[]) {
	
    if (argc == 1)
    {
		main_help();
        return (0);
    }
	
	
	
	
	homedir = getpwuid(getuid())->pw_dir;
	configfile = homedir + "/.phasex";
	load_config();
	

	int a;
	for (a = 2; a < argc; a++)
    {
		
        string str = argv[a];
		
        if (str.find("vcf=") != string::npos)
        {
            v_vcf = str.substr(4);
			continue;
        }
		else if (str.find("beagle=") != string::npos)
		{
			v_beagle = str.substr(7);
			continue;
		}
		else if (str.find("shapeit=") != string::npos)
		{
			v_shapeit = str.substr(8);
			continue;
		}

		else if (str.find("map=") != string::npos)
		{
			v_map = str.substr(4);
			continue;
		}
		else if (str.find("iterations=") != string::npos)
		{
			v_iteration = stoi(str.substr(11));
			if (v_iteration < 0) {v_iteration = 1;}
			continue;
		}
		else if (str.find("replicates=") != string::npos)
		{
			v_replicates = stoi(str.substr(11));
			if (v_replicates < 1) {v_replicates = 50;}
			continue;
		}
		
		else if (str.find("select=") != string::npos)
		{
			v_select = stof(str.substr(7));
			if (v_select > 1) {v_select = 0.7;}
			continue;
		}
		
		else if (str.find("threshold=") != string::npos)
		{
			v_threshold = stof(str.substr(10));
			if (stof(str.substr(10)) < 0.7) {v_threshold=0.7;}
			continue;
		}

		else if (str.find("minfix=") != string::npos)
		{
			min_fix = stoi(str.substr(7));
			if (stoi(str.substr(7)) < 1) {min_fix=1;}
			continue;
		}
		
		else if (str.find("input=") != string::npos)
		{
			v_input = str.substr(6);
			continue;
		}
		
		else if (str.find("output=") != string::npos)
		{
			v_output = str.substr(7);
			v_output = v_output + "/";
			continue;
		}
		
	
		else if (str.find("scheme=") != string::npos)
		{
			shapeit_scheme = str.substr(7);
			shapeit_scheme.erase(std::remove(shapeit_scheme.begin(), shapeit_scheme.end(), '"'), shapeit_scheme.end());
			shapeit_scheme.erase(std::remove(shapeit_scheme.begin(), shapeit_scheme.end(), '\''), shapeit_scheme.end());
			continue;
		}
		

		else if (str.find("shapeit_others=") != string::npos)
		{
			shapeit_other_parameters = str.substr(15);
			shapeit_other_parameters.erase(std::remove(shapeit_other_parameters.begin(), shapeit_other_parameters.end(), '"'), shapeit_other_parameters.end());
			shapeit_other_parameters.erase(std::remove(shapeit_other_parameters.begin(), shapeit_other_parameters.end(), '\''), shapeit_other_parameters.end());
			continue;
		}
		
		
		else if (str.find("threads=") != string::npos)
		{
			v_threads = str.substr(8);
			continue;
		}
		
		else if (str.find("--quiet") != string::npos)
		{
			v_quiet = 1;
			continue;
		}
		else if (str.find("--biallelic") != string::npos)
		{
			onlybiallelic = 1;
			continue;
		}
		else if (str.find("--debug") != string::npos)
		{
			debug = 1;
			continue;
		}
		/*



		else if (str.find("wh=") != string::npos)
		{
			v_ws = str.substr(3);
			continue;
		}

		else if (str.find("ws=") != string::npos)
		{
			v_ws = str.substr(3);
			continue;
		}
		
		else if (str.find("hc=") != string::npos)
		{
			v_hc = str.substr(3);
			continue;
		}
		
		else if (str.find("bin=") != string::npos)
		{
			v_phase = str.substr(4);
			continue;
		}
		
		else if (str.find("parameters=") != string::npos)
		{
			v_beagle_parameters = str.substr(11);
			v_phase_parameters = str.substr(11);
			continue;
		}

		else if (str.find("reduction=") != string::npos)
		{
			v_min_threshold = stof(str.substr(10));
			if (stof(str.substr(10)) < 0.7) {v_min_threshold=0.7;}
			continue;
		}




		else if (str.find("memory=") != string::npos)
		{
			v_memory = str.substr(7);
			continue;
		}


		else if (str.find("--unphase") != string::npos)
		{
			v_unphase = 1;
			continue;
		}

		
		else if (str.find("--clean") != string::npos)
		{
			v_cleanup = 1;
			continue;
		}

		else if (str.find("--debug") != string::npos)
		{
			v_debug = 1;
			continue;
		}

		else if (str.find("input=") != string::npos)
		{
			v_input = str.substr(6);
			continue;
		}

		else if (str.find("--nocombined") != string::npos)
		{
			v_nocombined = 1;
			continue;
		}
		*/
		else {
			v_message = "Unknown option: " + str;
			warnings.push_back(v_message);
		}
		
    }
	
	
	
	
    if (strcmp(argv[1],"phase-ps") == 0)
	{
		main_phase_ps();
	}

	else if (strcmp(argv[1],"recreate") == 0)
	{
		main_vcf();
	}

	else if (strcmp(argv[1],"hp-ps") == 0)
	{
		main_hp_ps_recode();
	}

	
	else
    {
		string str = argv[1];
		v_message = "Unknown command: " + str;
		warnings.push_back(v_message);
		main_help();
    }
	
	
	auto clock_end = std::chrono::steady_clock::now();
	auto diff = clock_end - clock_start;
	v_message = "Elapsed time: " + to_string(((std::chrono::duration <double, std::milli> (diff).count())/1000)) + " s";
	warnings.push_back(v_message);
	
	print_warnings();
	
    return 0;
}
