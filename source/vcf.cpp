
#include "external.h"
#include "functions.hpp"
#include <boost/algorithm/string.hpp>
#include <vector>
#include <map>


void vcf_help ()
{
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 0, Program_name + "::RECREATE", 1, 0);
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 2, "* Author  : " + Program_author, 1, 0);
    screen_message (screen_size, 2, "* Contact : " + Program_contact, 1, 0);
    screen_message (screen_size, 2, "* Version : " + Program_version, 1, 0);
    screen_message (screen_size, 2, "* Website : " + v_website, 1, 0);
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 2, "Options", 1, 0);
    screen_message (screen_size, 4, "input         Path to a phasex output folder", 1, 0);
    screen_message (screen_size, 4, "--quiet       quiet mode", 1, 0);
    screen_message (screen_size, 0, "", 1, 0);
    return;
    
}


void main_vcf ()
{
    if (v_input == "")  {warnings.push_back("You must indicate a phasex output folder!");vcf_help(); return;}
    if (! fileExists(v_input)) {warnings.push_back("This input phasex folder is not valid!");vcf_help(); return;}
    if (! fileExists(v_input + "/phasex.log")) {warnings.push_back("This input phasex folder is not valid!");vcf_help(); return;}

    string v_command = "grep Type " + v_input + "/phasex.log";
    string v_system_out = GetStdoutFromCommand(v_command);
    if (v_system_out == "") {warnings.push_back("This input phasex folder is not valid!");vcf_help(); return;}
    
    
    multiallelic = 0;
    size_t found=v_system_out.find("biallelic");
    if (found!=std::string::npos) {multiallelic = 0;}
    found=v_system_out.find("multi-allelic");
    if (found!=std::string::npos) {multiallelic = 1;}

    
    map <pair<string,int>,string> snpdata;
    sample_list.clear();
    int snpsize = 0;
    map <string, int> haplocount;
    int haplototal = 0;
    
    string input = "";
    if ((onlybiallelic == 0) && (multiallelic == 1)) {input = v_input + "/beagle/results.txt";}
    if (onlybiallelic == 1)  {input = v_input + "/shapeit/results.txt";}
    if (multiallelic == 0) {input = v_input + "/shapeit/results.txt";}

    
    
    ifstream inp (input);
    string line;
    while ( getline (inp,line))
    {
        if (line == "") {continue;}
        if (line.find("pass") != string::npos)
        {
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            string sample = data[0];
            sample_list.push_back(sample);
            string h1 = data[1];
            string h2 = data[2];
            
            haplocount[h1]++;
            haplocount[h2]++;
            haplototal++;
            haplototal++;
 
            vector <string> alleles_h1;
            vector <string> alleles_h2;
            boost::split(alleles_h1,h1,boost::is_any_of(" "));
            boost::split(alleles_h2,h2,boost::is_any_of(" "));
            
            for (int a = 0; a < alleles_h1.size(); a++)
            {
                pair<string, int> key;
                key = make_pair(sample,a);
                snpdata[key] = alleles_h1[a] + "|" + alleles_h2[a];
            }
        }
    }
    inp.close();
    
    
    
    string vcfinp = "";
    string test = v_input + "/source/original.vcf";
    if ((fileExists(test)) && (vcfinp == "")) {vcfinp = v_input + "/source/original.vcf";}
    //if (vcfinp == "") {cout << "Warning: this is not a phasex result folder!" << endl; return;}
    
    ifstream vcf (vcfinp);
    vector <string> snpinfo;
    while ( getline (vcf,line))
    {
        if (line.substr(0,1) == "#") {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        snpinfo.push_back(data[0] + ";" + data[1] + ";" + data[2] + ";" + data[3] + ";" + data[4]);
        snpsize++;
    }
    vcf.close();
    
    
    
    ofstream out;
    string outfile = v_input + "/results.vcf";
    out.open (outfile);

    out << "##fileformat=VCFv4.2\n";
    out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
    {out << "\t" << *samp;}
    out << "\n";

    
    for (int a = 0; a < snpsize; a++)
    {
        vector <string> data;
        boost::split(data,snpinfo[a],boost::is_any_of(";"));
        out << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t.\t.\t.\tGT";
        
        for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
        {
            pair<string, int> key;
            key = make_pair(*samp,a);
            out << "\t" << snpdata[key];
        }
        out << endl;
    }
    out.close();
    
    
    
    
    ofstream freq;
    outfile = v_input + "/results.freq";
    freq.open (outfile);
    
    freq << "Haplotype\tCount\tFrequency\n";
    
    for (auto &p : haplocount)
    {
        float frequency;
        frequency = (float)haplocount[p.first] / (float)haplototal;
        freq << p.first << "\t" << to_string(haplocount[p.first]) << "\t" << frequency << endl;
    }
    freq.close();

    
    
    outfile = v_input + "/sample_list.txt";
    freq.open (outfile);
    
    for(auto&& x: sample_list)
        freq << x << '\n';
    freq.close();

}
