#include <boost/algorithm/string.hpp>

#include "recode.hpp"
#include "functions.hpp"
#include "external.h"

void hp_ps_recode_help ()
{
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 0, Program_name + "::HP-PS", 1, 0);
    screen_message (screen_size, 0, "Recode ReadBackedPhasing VCF to be compatible with phasex", 1, 0);
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 2, "* Author  : " + Program_author, 1, 0);
    screen_message (screen_size, 2, "* Contact : " + Program_contact, 1, 0);
    screen_message (screen_size, 2, "* Version : " + Program_version, 1, 0);
    screen_message (screen_size, 2, "* Website : " + v_website, 1, 0);
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 2, "Options", 1, 0);
    screen_message (screen_size, 4, "vcf           VCF file phased using ReadBackedPhasing", 1, 0);
    screen_message (screen_size, 4, "output        The VCF file to be created", 1, 0);
    screen_message (screen_size, 4, "--quiet       quiet mode", 1, 0);
    screen_message (screen_size, 0, "", 1, 0);
    return;
}




void Recode_HP_PS (string input, string output)
{
    vector <string> samples;
    int sample_index = 0;
    ifstream inp (input);
    string line;
    
    ofstream outfile;
    string out = v_output;
    outfile.open (out);
    
    while ( getline (inp,line))
    {
        if (line.substr(0,2) == "##") {outfile << line << endl; continue;}
        if (line.substr(0,2) == "#C")
        {
            outfile << line << endl;
            boost::split(samples,line,boost::is_any_of("\t"));
            for(vector<string>::iterator iter = samples.begin();iter!=samples.end();++iter)
            {
                sample_index++;
                if (*iter == "FORMAT") {break;}
            }
            continue;
        }
        
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        vector <string> alt;
        boost::split(alt,data[4],boost::is_any_of(","));
        
        string format = data[sample_index-1];
        vector <string> itens;
        boost::split(itens,format,boost::is_any_of(":"));
        
        int GT = 0;
        int HP = 0;
        for(int a = 0; a < itens.size(); a++)
        {
            if (itens[a] == "GT") {GT = a;}
            if (itens[a] == "HP") {HP = a;}
        }
        
        for (int a = sample_index; a < data.size(); a++)
        {
            vector <string> values;
            boost::split(values,data[a],boost::is_any_of(":"));
            if (values[GT] == ".") {values[GT] = "./.";}
            string newdata;
            
            if (HP == 0)
            {
                newdata = values[GT] + ":.";
            }
            
            if ((HP != 0) && (values[HP] != "."))
            {
                
                string hp_value = values[HP];
                vector <string> sub;
                boost::split(sub,hp_value,boost::is_any_of(","));
                vector <string> vec;
                boost::split(vec,sub[0],boost::is_any_of("-"));
                string vec_data = vec[0];
                
                
                vector <string> alleles;
                boost::split(alleles,values[GT],boost::is_any_of("/"));
                if (alleles.size() == 1) {boost::split(alleles,values[GT],boost::is_any_of("|"));}

                string newgenotype;
                if (vec[1] == "1") {newgenotype = alleles[0] + "|" + alleles[1];}
                if (vec[1] == "2") {newgenotype = alleles[1] + "|" + alleles[0];}
                
                newdata = newgenotype + ":" + vec_data;
            }
           
            if ((HP != 0) && (values[HP] == "."))
            {
                newdata = values[GT] + ":.";
            }
            
            data[a] = newdata;
        }
        
        
   
            for (int a = 0; a < sample_index - 1; a++)
            {
                outfile << data[a] << "\t";
            }
            outfile << "GT:PS";
        
            for (int a = sample_index; a < data.size(); a++)
            {
                outfile << "\t" << data[a];
            }
            outfile << endl;
     
    }
    outfile.close();
    inp.close();
}








void main_hp_ps_recode ()
{
    if (v_vcf == "")  {warnings.push_back("You must indicate a VCF file!");hp_ps_recode_help(); return;}
    if (! fileExists(v_vcf)) {warnings.push_back("You must indicate a VCF file!");hp_ps_recode_help(); return;}
    if (v_output == "")  {warnings.push_back("You must indicate the VCF file to be created!");hp_ps_recode_help(); return;}
    
    
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 0, Program_name + "::HP-PS", 1, v_quiet);
    
    screen_message (screen_size, 2, "Recoding file ...", 2, v_quiet);

    v_output = v_output.substr(0,v_output.size() - 1);
    
    
    Recode_HP_PS (v_vcf, v_output);
    
    screen_message (screen_size, 2, "Recoding file ... done", 1, v_quiet);
    screen_message (screen_size, 2, "", 1, v_quiet);

    return;
}




