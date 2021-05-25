#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <vector>
#include "ThreadPool.h"
#include <mutex>

#include "functions.hpp"
#include "external.h"


std::mutex mtx_read_result;

using namespace std;


void print_debug (int debug, string msg)
{
    if (debug == 0) {return;}
    cout << "   Debug >> " << msg << endl;
}



void print_warnings ()
{
    if (warnings.size() > 0)
    {
        v_message = "Warning:  " + warnings[0];
        screen_message (screen_size, 2, v_message, 1, v_quiet);
    }
    if (warnings.size() > 1)
    {
        for (int a = 1; a < warnings.size(); a++) {
            v_message = "          " + warnings[a];
            screen_message (screen_size, 2, v_message, 1, v_quiet);
        }
    }
    screen_message (screen_size, 0, "", 1, v_quiet);
}


void screen_message (int size, int left, string message, int enter, int quiet)
{
    if(quiet == 1) {return;
        
    }
    if ((message.length()+left) > size) {message = message.substr(0,(size-left-1));}
    cout << "\r";
    for (int a = 0; a < left; a++){cout << " ";}
    cout << message;
    for (int a = 0; a < (((size-left)-message.length())); a++){cout << " ";}
    std::cout.flush();
    if (enter == 1) {cout << endl;}
    return;
}


bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

string findfilepath (string v_file)
{
    string v_filename = v_file.substr(0,(v_file.find_last_of("/"))+1);
    return (v_filename);
}

string GetStdoutFromCommand(string cmd) {
    
    string data;
    FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    cmd.append(" 2>&1");
    
    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
        pclose(stream);
    }
    return data;
}

void printmessage (string v_text, int v_quiet, int v_level, int v_enter)
{
    if (v_quiet == 1) {return;}
    
    for (int v_temp = 0; v_temp < v_level; v_temp++)
    {
        cout << "\t";
    }
    cout << v_text;
    if (v_enter == 1) {cout << endl;}
    return;
}





bool ends_with(const std::string &filename, const std::string &ext)
{
    return ext.length() <= filename.length() &&
    std::equal(ext.rbegin(), ext.rend(), filename.rbegin());
}





void PrepareVcfPS (string file)
{
    vector <string> samples;
    int sample_index = 0;
    ifstream inp (file);
    string line;
    
    ofstream biallele;
    string outbi = v_output + "/source/biallelic.vcf";
    biallele.open (outbi);

    ofstream psallele;
    string outps = v_output + "/source/original-ps.vcf";
    psallele.open (outps);

    ofstream multiallele;
    string outmulti = v_output + "/source/multiallelic.vcf";
    multiallele.open (outmulti);

    while ( getline (inp,line))
    {
        if (line.substr(0,2) == "##") {biallele << line << endl; multiallele << line << endl; psallele << line << endl; continue;}
        if (line.substr(0,2) == "#C")
        {
            biallele << line << endl; multiallele << line << endl; psallele << line << endl;
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
        int PS = 0;
        for(int a = 0; a < itens.size(); a++)
        {
            if (itens[a] == "GT") {GT = a;}
            if (itens[a] == "PS") {PS = a;}
        }
        
        for (int a = sample_index; a < data.size(); a++)
        {
            vector <string> values;
            boost::split(values,data[a],boost::is_any_of(":"));
            if (values[GT] == ".") {values[GT] = "./.";}
            string newdata;
            if (PS == 0) {newdata = values[GT] + ":.";}
            if (PS != 0) {newdata = values[GT] + ":" + values[PS];}
            data[a] = newdata;
        }
        
        
        if (alt.size() == 1) {
            for (int a = 0; a < sample_index - 1; a++)
            {
                biallele << data[a] << "\t";
                psallele << data[a] << "\t";
            }
            biallele << "GT:PS";
            psallele << "GT:PS";
            
            for (int a = sample_index; a < data.size(); a++)
            {
                biallele << "\t" << data[a];
                psallele << "\t" << data[a];
            }
            biallele << endl;
            psallele << endl;
        }
        

        if (alt.size() > 1) {
            multiallelic = 1;
            for (int a = 0; a < sample_index - 1; a++)
            {
                multiallele << data[a] << "\t";
                psallele << data[a] << "\t";
            }
            multiallele << "GT:PS";
            psallele << "GT:PS";
            
            for (int a = sample_index; a < data.size(); a++)
            {
                multiallele << "\t" << data[a];
                psallele << "\t" << data[a];
            }
            multiallele << endl;
            psallele << endl;
        }
        
    }
    multiallele.close();
    biallele.close();
    psallele.close();
    inp.close();
    if (multiallelic == 0){onlybiallelic = 1;}
}






string CheckVCF (string file)
{
    snp.clear();
    int sample_index = 0;
    vector <string> samples;
    map <string,int> regions;
    int sample_size = 0;
    int snp_bi = 0;
    int snp_mm = 0;
    int gtps = 1;
    int HP = -1;

    
    ifstream inp (file);
    string line;
    while ( getline (inp,line))
    {
        if (line.substr(0,2) == "##") {continue;}
        if (line.substr(0,2) == "#C")
        {
            boost::split(samples,line,boost::is_any_of("\t"));
            for(vector<string>::iterator iter = samples.begin();iter!=samples.end();++iter)
            {
                sample_index++;
                if (*iter == "FORMAT") {break;}
            }
            for (int a = sample_index; a < samples.size(); a++)
            {
                sample_size++;
            }
            continue;
        }
        
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        
        
        string format = data[sample_index-1];
        vector <string> itens;
        boost::split(itens,format,boost::is_any_of(":"));
        
        int GT = -1;
        int PS = -1;
        for(int a = 0; a < itens.size(); a++)
        {
            if (itens[a] == "GT") {GT = a;}
            if (itens[a] == "PS") {PS = a;}
            if (itens[a] == "HP") {HP = a;}
        }
        
        if ((GT == -1) || (PS == -1)){gtps = 0; break;}
        
        regions[data[0]] = 1;
        std::size_t found = data[4].find(",");
        if (found!=std::string::npos) {snp_mm++;}
        else {snp_bi++;}
    }
    inp.close();

    if (gtps == 0)
    {
        warnings.push_back("Some variants do not present the GT or PS fields.");
        if (HP != -1) {
            warnings.push_back("This seens to be a ReadBackedPhasing VCF.");
            warnings.push_back("Consider recode this with 'phasex hp-ps'.");
        }
        return "";
    }
    
    int count_regions = 0;
    for (auto const& pair: regions)
    {
        count_regions++;
        v_region = pair.first;
    }
    if (count_regions > 1) {v_region = "";warnings.push_back("More than one chromosome. This is no compatible with phasex.");return "";}
    
    string res = to_string(sample_size) + " samples, " + to_string(snp_bi) + " biallelic variants, " + to_string(snp_mm) + " multi-allelic variants, region " + v_region;
    return res;
}




void LoadingGenotypes (string file)
{
    db.clear();
    snp.clear();
    int sample_index = 0;
    sample_list.clear();
    snp_list.clear();
    vector <string> samples;
    
    ifstream inp (file);
    string line;
    while ( getline (inp,line))
    {
        if (line.substr(0,2) == "##") {continue;}
        if (line.substr(0,2) == "#C")
        {
            boost::split(samples,line,boost::is_any_of("\t"));
            for(vector<string>::iterator iter = samples.begin();iter!=samples.end();++iter)
            {
                sample_index++;
                if (*iter == "FORMAT") {break;}
            }
            for (int a = sample_index; a < samples.size(); a++)
            {
                sample_list.push_back(samples[a]);
            }
            
            continue;
        }
        
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        
        string keysnp = data[1] + "," + data[3] + "," + data[4];
        snp_list.push_back(keysnp);
        snp[keysnp] = data[0] + ";" + data[1] + ";" + data[2] + ";" + data[3] + ";" + data[4];
        
        for (int a = sample_index; a < data.size(); a++)
        {
            pair<string, string> key;
            key = make_pair(samples[a],keysnp);
            db[key] = data[a];
        }
    }
    inp.close();
}


void LoadingGenotypesMulti (string file)
{
    int sample_index = 0;
    vector <string> samples;
    
    ifstream inp (file);
    string line;
    while ( getline (inp,line))
    {
        if (line.substr(0,2) == "##") {continue;}
        if (line.substr(0,2) == "#C")
        {
            boost::split(samples,line,boost::is_any_of("\t"));
            for(vector<string>::iterator iter = samples.begin();iter!=samples.end();++iter)
            {
                sample_index++;
                if (*iter == "FORMAT") {break;}
            }
            for (int a = sample_index; a < samples.size(); a++)
            {
                //sample_list.push_back(samples[a]);
            }
            
            continue;
        }
        
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        
        string keysnp = data[1] + "," + data[3] + "," + data[4];

        snp_list.push_back(keysnp);
        snp[keysnp] = data[0] + ";" + data[1] + ";" + data[2] + ";" + data[3] + ";" + data[4];
        
        for (int a = sample_index; a < data.size(); a++)
        {
            pair<string, string> key;
            key = make_pair(samples[a],keysnp);
            db[key] = data[a];
        }
    }
    inp.close();
    
    sort(snp_list.begin(), snp_list.end());
    sort(sample_list.begin(), sample_list.end());

}





void ReadResults (string tag, int replicates, int iteration)
{
    final_haplotypes.clear();
    
    ThreadPool pool(stoi(v_threads));
    std::vector< std::future<int> > results;

    
    for (int replica = 1; replica <= replicates; replica++)
    {
        results.emplace_back(
            pool.enqueue([replica, replicates, iteration, tag] {
                map <string, string> haplo_sample_h1;
                map <string, string> haplo_sample_h2;
                vector <string> samples;
                string file = v_output + "/" + tag + "/iteration_" + to_string(iteration) + "/phased_" + to_string(replica) + ".vcf";
                ifstream inp (file);
                string line;
            
                while ( getline (inp,line))
                {
                    if (line.substr(0,2) == "##") {continue;}
                    if (line.substr(0,2) == "#C")
                    {
                        boost::split(samples,line,boost::is_any_of("\t"));
                        continue;
                    }
                    vector <string> data;
                    boost::split(data,line,boost::is_any_of("\t"));
                    for (int b = 9; b < data.size();b++)
                    {
                        vector <string> alleles;
                        boost::split(alleles,data[b],boost::is_any_of("|"));
                        haplo_sample_h1[samples[b]] = haplo_sample_h1[samples[b]] + " " + alleles[0];
                        haplo_sample_h2[samples[b]] = haplo_sample_h2[samples[b]] + " " + alleles[1];
                    }
                }
                inp.close();
            
                for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
                {
                    string finalpair = haplo_sample_h1[*samp].substr(1) + "\t" + haplo_sample_h2[*samp].substr(1);
                    string finalpair_rev = haplo_sample_h2[*samp].substr(1) + "\t" + haplo_sample_h1[*samp].substr(1);
                    
                    if (used_pairs.count(finalpair) > 0)
                    {
                        pair<string, string> key;
                        key = make_pair(*samp,finalpair);
                        mtx_read_result.lock();
                        final_haplotypes[key]++;
                        mtx_read_result.unlock();
                        continue;
                    }
                    if (used_pairs.count(finalpair_rev) > 0)
                    {
                        pair<string, string> key;
                        key = make_pair(*samp,finalpair_rev);
                        mtx_read_result.lock();
                        final_haplotypes[key]++;
                        mtx_read_result.unlock();
                        continue;
                    }
                    pair<string, string> key;
                    key = make_pair(*samp,finalpair);
                    mtx_read_result.lock();
                    final_haplotypes[key]++;
                    used_pairs[finalpair] = 1;
                    mtx_read_result.unlock();
                }
            
            
                return replica*replica;
            })
        );

    } //end replicates
    for(auto && result: results){result.get();} // waiting for all threads
    

   if ((tag == "shapeit") || (tag == "beagle"))
    {
        final_haplotypes_vcf.clear();
        final_haplotypes_vcf = final_haplotypes;
    }
}



void LoadingGenotypesFixed ()
{

    for (auto &p : fixed_haplotypes)
    {
        print_debug (debug, "LoadingGenotypeFixed: sample " + p.first);
        string samp = p.first;
        string haplo = fixed_haplotypes[samp];
        vector <string> h;
        boost::split(h,haplo,boost::is_any_of("\t"));
        vector <string> nucA;
        vector <string> nucB;
        boost::split(nucA,h[0],boost::is_any_of(" "));
        boost::split(nucB,h[1],boost::is_any_of(" "));
        int vec = 0;
        string snp_first = "";
        for(vector<string>::iterator it = snp_list.begin();it!=snp_list.end();++it)
        {
            if (snp_first == "")
            {
                vector <string> keysnp_split;
                boost::split(keysnp_split,*it,boost::is_any_of(","));
                snp_first = keysnp_split[0];
            }
            string gen = nucA[vec] + "|" + nucB[vec] + ":" + snp_first;
            pair<string, string> key;
            key = make_pair(samp,*it);
            db[key] = gen;
            vec++;
        }
        
    }
}




void PrintVCFPS(string tag, int iter)
{
    ofstream myfile;
    string outfile;
    if (tag == "shapeit") {outfile = v_output + "/" + tag + "/iteration_" + to_string(iter) + "/phased_final.vcf";}
    if ((tag == "beagle") && (iter == 0)) {outfile = v_output + "/" + tag + "/source/source.vcf";}
    if ((tag == "beagle") && (iter > 0)) {outfile = v_output + "/" + tag + "/iteration_" + to_string(iter) + "/phased_final.vcf";}
    myfile.open (outfile);
    
    myfile << "##fileformat=VCFv4.2\n";
    myfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    myfile << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n";

    myfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
    {myfile << "\t" << *samp;}
    myfile << endl;
    
    for(vector<string>::iterator it = snp_list.begin();it!=snp_list.end();++it)
    {
        vector <string> data;
        boost::split(data,snp[*it],boost::is_any_of(";"));
        myfile << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t";
        myfile <<".\t.\t.\tGT:PS";
        
        string position = *it;
        
        for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
        {
            string sample = *samp;
            
            pair<string, string> key;
            key = make_pair(*samp,position);
            string snp_data = db[key];
            
            vector <string> alleles;
            boost::split(alleles,snp_data,boost::is_any_of("/"));
            int phased = 0;
            if (alleles.size() == 1) {boost::split(alleles,snp_data,boost::is_any_of("|"));phased = 1;}
            
            
            //if ((find(list_next_round.begin(), list_next_round.end(), *samp) != list_next_round.end()) || (phased == 1))
            if ((find(list_next_round.begin(), list_next_round.end(), *samp) != list_next_round.end()))
            {
                if (phased == 1) {
                    myfile << "\t" << alleles[0] << "|" << alleles[1];
                    continue;
                }
                if (phased == 0) {
                    myfile << "\t" << alleles[0] << "/" << alleles[1];
                    continue;
                }
            }
            
            myfile << "\t" << alleles[0] << "/" << alleles[1];
            continue;
        }
        myfile << "\n";
    }
    myfile.close();
    
    
    if (tag == "shapeit") {
        outfile = v_output + "/" + tag + "/iteration_" + to_string(iter) + "/phased_final.vcf";
        string v_command = "bgzip -c  " + outfile + " > " + outfile + ".gz";
        string v_system_out = GetStdoutFromCommand(v_command);
        v_command = "tabix -p vcf " + outfile + ".gz";
        v_system_out = GetStdoutFromCommand(v_command);
    }
}



void PrintVCFPSnew(string tag, int iter)
{
    ofstream myfile;
    string outfile;
    if (tag == "shapeit") {outfile = v_output + "/" + tag + "/iteration_" + to_string(iter) + "/phased_final.vcf";}
    if ((tag == "beagle") && (iter == 0)) {outfile = v_output + "/" + tag + "/source/source.vcf";}
    if ((tag == "beagle") && (iter > 0)) {outfile = v_output + "/" + tag + "/iteration_" + to_string(iter) + "/phased_final.vcf";}
    myfile.open (outfile);
    
    myfile << "##fileformat=VCFv4.2\n";
    myfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    myfile << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n";

    myfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
    {myfile << "\t" << *samp;}
    myfile << endl;
    
    map <string,string> vcf_data_thread;
    ThreadPool pool(stoi(v_threads));
    std::vector< std::future<int> > results;

    int count = 0;
    for(vector<string>::iterator it = snp_list.begin();it!=snp_list.end();++it)
    {
        string snploop = *it;
        results.emplace_back(
            pool.enqueue([count,snploop, &vcf_data_thread] {
            
            string line = "";
            vector <string> data;
            boost::split(data,snp[snploop],boost::is_any_of(";"));
            line =  data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + data[3] + "\t" + data[4] + "\t";
            line.append(".\t.\t.\tGT:PS");
        
            string position = snploop;
        
            for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
            {
                string sample = *samp;
            
                pair<string, string> key;
                key = make_pair(*samp,position);
                string snp_data = db[key];
            
                vector <string> alleles;
                boost::split(alleles,snp_data,boost::is_any_of("/"));
                int phased = 0;
                if (alleles.size() == 1) {boost::split(alleles,snp_data,boost::is_any_of("|"));phased = 1;}
            
            
                if ((find(list_next_round.begin(), list_next_round.end(), *samp) != list_next_round.end()))
                {
                    if (phased == 1) {
                        line.append("\t" + alleles[0] + "|" + alleles[1]);
                        continue;
                    }
                    if (phased == 0) {
                        line.append("\t" + alleles[0] + "/" + alleles[1]);
                        continue;
                    }
                }
            
                line.append("\t" + alleles[0] + "/" + alleles[1]);
                continue;
        }
        mtx_read_result.lock();
        vcf_data_thread[snploop] = line;
        mtx_read_result.unlock();
        return 1;
        })
        );
    }
    for(auto && result: results){result.get();} // waiting for all threads
    
    for(vector<string>::iterator it = snp_list.begin();it!=snp_list.end();++it)
    {
        myfile << vcf_data_thread[*it] << endl;
    }
    myfile.close();
    
    
    if (tag == "shapeit") {
        outfile = v_output + "/" + tag + "/iteration_" + to_string(iter) + "/phased_final.vcf";
        string v_command = "bgzip -c  " + outfile + " > " + outfile + ".gz";
        string v_system_out = GetStdoutFromCommand(v_command);
        v_command = "tabix -p vcf " + outfile + ".gz";
        v_system_out = GetStdoutFromCommand(v_command);
    }
}



void PrintVCFBeagle()
{
    
    
    ofstream myfile;
    string outfile = v_output + "/beagle/source/source.vcf";
    myfile.open (outfile);
    myfile << "##fileformat=VCFv4.2\n";
    myfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    myfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp) {myfile << "\t" << *samp;}
    myfile << endl;
    
    for(vector<string>::iterator it = snp_list.begin();it!=snp_list.end();++it)
    {
        vector <string> data;
        boost::split(data,snp[*it],boost::is_any_of(";"));
        myfile << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t";
        myfile <<".\t.\t.\tGT";
        
        string position = *it;
        
        for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
        {
            string sample = *samp;
            
            pair<string, string> key;
            key = make_pair(*samp,position);
            string snp_data = db[key];
            
            vector <string> data;
            boost::split(data,snp_data,boost::is_any_of(":"));
            vector <string> alleles;
            boost::split(alleles,data[0],boost::is_any_of("/"));
            int phased = 0;
            if (alleles.size() == 1) {boost::split(alleles,data[0],boost::is_any_of("|"));phased = 1;}
            
            
            if ((find(list_next_round.begin(), list_next_round.end(), *samp) != list_next_round.end()))
            {
                if (phased == 1) {
                    myfile << "\t" << alleles[0] << "|" << alleles[1];
                    continue;
                }
                if (phased == 0) {
                    myfile << "\t" << alleles[0] << "/" << alleles[1];
                    continue;
                }
            }
            
            myfile << "\t" << alleles[0] << "/" << alleles[1];
            continue;
        }
        myfile << "\n";
    }
    myfile.close();
     
}


/*


void CreateOriginalKnown(string file, string tag, string out)
{
    vector<string> samples;
    int sample_index = 0;
    map <pair<string,string>, string> data;
    map <pair<string,string>, string> known;
    map <pair<string,string>, int> hp_db;
    vector <string> positions;
    vector <string> sample_list;
    
    
    string line;
    ifstream knownfile (file);
    while ( getline (knownfile,line) )
    {
        if (line.substr(0,2) == "##") {continue;}
        if (line.substr(0,2) == "#C")
        {
            boost::split(samples,line,boost::is_any_of("\t"));
            for(vector<string>::iterator it = samples.begin();it!=samples.end();++it)
            {
                sample_index++;
                if (*it == "FORMAT"){break;}
            }
            continue;
        }
        
        vector <string> linedata;
        boost::split(linedata,line,boost::is_any_of("\t"));
        positions.push_back(linedata[1]);
        int hp = 0;
        int hp_found = 0;
        vector <string> hps;
        boost::split(hps,linedata[(sample_index-1)],boost::is_any_of(":"));
        
        for(vector<string>::iterator it = hps.begin();it!=hps.end();++it)
        {
            if (*it == "HP"){hp_found = 1;break;}
            hp++;
        }
        if (hp == hps.size()) {hp = 1000;}
        
        for (int a = sample_index; a < linedata.size(); a++)
        {
            sample_list.push_back(samples[a]);
            pair<string, string> key;
            key = make_pair(samples[a],linedata[1]);
            data[key] = linedata[a];
            hp_db[key] = hp;
        }
    }
    knownfile.close();
    
    
    
    for (int a = sample_index; a < samples.size(); a++)
    {
        
        int first = 0;
        map <string, int> vectors;
        
        for(vector<string>::iterator pos = positions.begin();pos!=positions.end();++pos)
        {
            string real_position = *pos;
            vector <string> snp_data;
            pair<string, string> key;
            key = make_pair(samples[a],real_position);
            boost::split(snp_data,data[key],boost::is_any_of(":"));
            vector <string> alleles;
            boost::split(alleles,snp_data[0],boost::is_any_of("/"));
            if (alleles.size() == 0) {boost::split(snp_data,snp_data[0],boost::is_any_of("|"));}
            
            if ((alleles[0] == ".") || (alleles[1] == ".")) {known[key] = "*";continue;}
//            if (alleles[0] == alleles[1]) {known[key] = "*";continue;}
            
            int hp_db_value = hp_db[key];
            if (hp_db_value == 1000) {known[key] = "*";continue;}
            if (snp_data[hp_db_value] == "."){known[key] = "*";continue;}
            
            if (hp_db_value != 1000)
            {
                
                vector <string> phases;
                
                boost::split(phases,snp_data[hp_db[key]],boost::is_any_of(","));
                
                vector <string> innervector;
                boost::split(innervector,phases[0],boost::is_any_of("-"));
                string refered_to = innervector[0];
                int phase = stoi(innervector[1]);
                
                
                
                if (refered_to == real_position)
                {
                    vectors[real_position]++;
                    if (first == 1) {
                        known[key] = "|" + to_string(phase - 1);
                    }
                    
                    if (first == 0) {
                        known[key]= to_string(phase - 1);
                        first = 1;
                    }
                    continue;
                }
                
                if (refered_to != real_position)
                {
                    first = 1;
                    vectors[refered_to]++;
                    known[key]= to_string(phase - 1);
                    continue;
                }
                
                continue;
                
            }
            
        }
        
        for(std::map<string,int>::iterator it = vectors.begin(); it != vectors.end(); ++it)
        {
            string pos = it->first;
            if (vectors[pos] == 1)
            {
                pair<string, string> key;
                key = make_pair(samples[a],pos);
                known[key] = "*";
            }
        }
        
    }
    
    
    
    ofstream myfile;
    string outfile = out + "/" + tag + ".known";
    myfile.open (outfile);
    
    for (int a = sample_index; a < samples.size(); a++)
    {
        myfile << samples[a] + ":";
        string line = "";
        for(vector<string>::iterator pos = positions.begin();pos!=positions.end();++pos)
        {
            string position = *pos;
            pair<string, string> key;
            key = make_pair(samples[a],position);
            line = line + known[key];
        }
        
        vector <string> k;
        boost::split(k,line,boost::is_any_of("|"));
        line = "";
        for(vector<string>::iterator ksub = k.begin();ksub!=k.end();++ksub)
        {
            string sub = *ksub;
            size_t found = sub.find("0");
            if (found != string::npos)
            {
                line = line + sub + "|";
                continue;
            }
            found = sub.find("1");
            if (found != string::npos)
            {
                line = line + sub + "|";
                continue;
            }
            
            line = line + sub;
        }
        
        if (line.substr((line.length()-1),1) == "|") {
            line = line.substr(0,line.length()-1);
        }
        myfile << line << endl;
    }
    myfile.close();
    return;
}

void LoadHaplotypeBlocks(string file, string item)
{
    string line;
    ifstream knownfile (file);
    
    if (item == "hc") {
        haplotype_block_hc.clear();
    }
    
    while ( getline (knownfile,line) )
    {
        vector <string> data;
        boost::split(data,line,boost::is_any_of(";"));
        
        if (item == "hc") {
            haplotype_block_hc[data[0]] = haplotype_block_hc[data[0]] + ";" + data[1];
        }
        
    }
    knownfile.close();
}

 
 void CreateOriginalKnownHC (string file, string tag, string out)
 {
     vector<string> samples;
     int sample_index = 0;
     map <pair<string,string>, string> data;
     map <pair<string,string>, string> known;
     map <pair<string,string>, int> hp_db;
     vector <string> positions;
     vector <string> sample_list;
 
 
 
     string line;
     ifstream knownfile (file);
     while ( getline (knownfile,line) )
     {
        if (line.substr(0,2) == "##") {continue;}
        if (line.substr(0,2) == "#C")
        {
            boost::split(samples,line,boost::is_any_of("\t"));
            for(vector<string>::iterator it = samples.begin();it!=samples.end();++it)
            {
                sample_index++;
                if (*it == "FORMAT"){break;}
            }
            continue;
        }
 
         vector <string> linedata;
         boost::split(linedata,line,boost::is_any_of("\t"));
         positions.push_back(linedata[1]);
         int hp = 0;
         int hp_found = 0;
         vector <string> hps;
         boost::split(hps,linedata[(sample_index-1)],boost::is_any_of(":"));
 
         for(vector<string>::iterator it = hps.begin();it!=hps.end();++it)
         {
             if (*it == "PID"){hp_found = 1;break;}
             hp++;
         }
         if (hp == hps.size()) {hp = 1000;}
 
         for (int a = sample_index; a < linedata.size(); a++)
         {
             sample_list.push_back(samples[a]);
             pair<string, string> key;
             key = make_pair(samples[a],linedata[1]);
             data[key] = linedata[a];
             hp_db[key] = hp;
         }
     }
     knownfile.close();
 
 
     ofstream myfile;
     string outfile = out + "/" + tag + ".known";
     myfile.open (outfile);
     
 
     for (int a = sample_index; a < samples.size(); a++)
     {
 
         int first = 0;
         map <string, string> vectors;
 
         for(vector<string>::iterator pos = positions.begin();pos!=positions.end();++pos)
         {
             string real_position = *pos;
             vector <string> snp_data;
             pair<string, string> key;
             key = make_pair(samples[a],real_position);
             boost::split(snp_data,data[key],boost::is_any_of(":"));
             vector <string> alleles;
             boost::split(alleles,snp_data[0],boost::is_any_of("/"));
             if (alleles.size() == 1) {boost::split(alleles,snp_data[0],boost::is_any_of("|"));}
             if ((alleles[0] == ".") || (alleles[1] == ".")) {continue;}
 
             int hp_db_value = hp_db[key];
             if (hp_db_value == 1000) {continue;}
             if (snp_data[hp_db_value] == "."){continue;}
 
             if (hp_db_value != 1000)
             {
                 vector <string> phases;
                 boost::split(phases,snp_data[hp_db[key]],boost::is_any_of("_"));
                 string refered_to = phases[0];
                 vectors[refered_to] = vectors[refered_to] + real_position + ":" + alleles[0] + ",";
 
                 continue;
             }
         }
 
         for(std::map<string,string>::iterator it = vectors.begin(); it != vectors.end(); ++it)
         {
             string pos = it->first;
             myfile << samples[a] << ";" + vectors[pos] << endl;
         }
 
     }
 
 
     myfile.close();
     return;
 }




void CombineBlocks (string file1, string file2, string out)
{
    ifstream known1 (file1);
    ifstream known2 (file2);
    ofstream myfile;
    string outfile = out + "/combined.known";
    myfile.open (outfile);
    
    string line1 = "";
    string line2 = "";
    while ( getline (known1,line1))
    {
        getline(known2, line2);
        vector <string> data;
        boost::split(data,line1,boost::is_any_of(":"));
        string sample = data[0];
        line1 = data[1];
        boost::split(data,line2,boost::is_any_of(":"));
        line2 = data[1];
        string combined = "";
        
        int max_size = max (line1.length(), line2.length());
        for (int a = 0; a < max_size; a++)
        {
            if ((line1.substr(a,1) == "|") && (line2.substr(a,1) != "|")) {line2.insert(a,"|");}
            if ((line2.substr(a,1) == "|") && (line1.substr(a,1) != "|")) {line1.insert(a,"|");}
        }
        
        
        vector <string> k1;
        vector <string> k2;
        boost::split(k1,line1,boost::is_any_of("|"));
        boost::split(k2,line2,boost::is_any_of("|"));
        
        myfile << sample << ":";
        
        
        for (int a = 0; a < k1.size(); a++)
        {
            string k1i = k1[a];
            string k2i = k2[a];
            
            int divergence_A = 0;
            for (int b = 0; b < k1i.length(); b++)
            {
                if ((k1i.substr(b,1) == "*") || (k2i.substr(b,1) == "*")){continue;}
                if (k1i.substr(b,1) != k2i.substr(b,1)){divergence_A++;continue;}
            }
            
            int divergence_B = 0;
            for (int b = 0; b < k1i.length(); b++)
            {
                if ((k1i.substr(b,1) == "*") || (k2i.substr(b,1) == "*")){continue;}
                if (k1i.substr(b,1) == k2i.substr(b,1)){divergence_B++;continue;}
            }
            
            if (divergence_B < divergence_A)
            {
                replace( k2i.begin(), k2i.end(), '0', 'b');
                replace( k2i.begin(), k2i.end(), '1', 'a');
                replace( k2i.begin(), k2i.end(), 'b', '1');
                replace( k2i.begin(), k2i.end(), 'a', '0');
            }
            
            
            for (int b = 0; b < k1i.length(); b++)
            {
                if ((k1i.substr(b,1) == "*") || (k2i.substr(b,1) == "*")){combined = combined + "*";continue;}
                if (k1i.substr(b,1) != k2i.substr(b,1)){combined = combined + "*";continue;}
                if (k1i.substr(b,1) == k2i.substr(b,1)){combined =  combined + k1i.substr(b,1);continue;}
            }
            combined = combined + "|";
        }
        
        boost::split(k1,combined,boost::is_any_of("|"));
        combined = "";
        
        for(vector<string>::iterator iter = k1.begin();iter!=k1.end();++iter)
        {
            string data = *iter;
            size_t found = data.find("0");
            if (found != string::npos)
            {
                combined = combined + data + "|";
                continue;
            }
            found = data.find("1");
            if (found != string::npos)
            {
                combined = combined + data + "|";
                continue;
            }
            combined = combined + data;
        }
        
        
        if (combined.substr((combined.length()-1),1) == "|") {
            combined = combined.substr(0,combined.length()-1);
        }
        myfile << combined << endl;
    }
    
    known1.close();
    known2.close();
    myfile.close();
}




int CompareHaploBlock(string sample, string haplo, string item){

    string data;
    if (item == "hc") {data = haplotype_block_hc[sample];}
    vector <string> blocks;
    boost::split(blocks,data,boost::is_any_of(";"));
    vector <string> haplos;
    boost::split(haplos,haplo,boost::is_any_of("\t"));
    vector <string> haplo1_data;
    boost::split(haplo1_data,haplos[0],boost::is_any_of(" "));
    vector <string> haplo2_data;
    boost::split(haplo2_data,haplos[1],boost::is_any_of(" "));

    map <string,int> conversion;
    for (int a = 0; a < snp_list.size(); a++)
    {
        conversion[snp_list[a]] = a;
    }
    
    int res = 0;
    int count_h1 = 0;
    int count_h2 = 0;
    int count = 0;
    
    for(vector<string>::iterator iter = blocks.begin();iter!=blocks.end();++iter)
    {
        if (*iter == "") {continue;}
        
        vector <string> subdata;
        boost::split(subdata,*iter,boost::is_any_of(","));
        if (subdata.size() == 2) {continue;}
        
        for(vector<string>::iterator sub = subdata.begin();sub!=subdata.end();++sub)
        {
            if (*sub == "") {continue;}
            
            vector <string> current;
            boost::split(current,*sub,boost::is_any_of(":"));
            string pos = current[0];
            string allele = current[1];
        
            if (haplo1_data[conversion[pos]] == allele) {count_h1++;}
            if (haplo2_data[conversion[pos]] == allele) {count_h2++;}
            count++;
        }
    }
    if (count_h1 == count) {res = 1;}
    if (count_h2 == count) {res = 1;}
    
    return res;
}



void LoadingBlocks (string file)
{
    
    knowndb.clear();
    knownblocks.clear();
    
    ifstream inp (file);
    string line;
    while ( getline (inp,line))
    {
        vector <string> data;
        boost::split(data,line,boost::is_any_of(":"));
        string k = data[1];
        k.erase(std::remove(k.begin(), k.end(), '|'), k.end());
        
        for (int a = 0; a < k.length(); a++)
        {
            pair<string, int> key;
            key = make_pair(data[0],a);
            knowndb[key] = k.substr(a,1);
        }

        size_t found = data[1].find("|");
        if (found == string::npos)
        {
            knownblocks[data[0]] = knownblocks[data[0]] + ";0," + to_string(k.length());
        }
        
        if (found != string::npos)
        {
            vector <string> block;
            boost::split(block,data[1],boost::is_any_of("|"));
            int vec = 0;
            for(vector<string>::iterator iter = block.begin();iter!=block.end();++iter)
            {
                string current = *iter;
                knownblocks[data[0]] =  knownblocks[data[0]] + ";" + to_string(vec) + "," + to_string(current.length());
                vec = vec + current.length();
            }
        }
    }
    inp.close();
}






void PrintingVCFs (int replicates, string tag, string out)
{
    for (int a = 1; a <= replicates; a++)
    {
        if (v_quiet == 0) {cout << "      > Writing VCF, replicate number " + to_string(a) + "\r" << flush;}
        
        randblock.clear();
        
        for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
        {
            vector <string> known;
            string knowndata = knownblocks[*samp].substr(1);
            boost::split(known,knowndata,boost::is_any_of(";"));
            if (known.size() == 1) {randblock[*samp] = 0;continue;}
            if (known.size() > 1)
            {
                int random = rand() % known.size();
                randblock[*samp] = random;
                continue;
            }
            
        }

        ofstream myfile;
        string outfile = out + "/replicate_" + to_string(a) + ".vcf";
        myfile.open (outfile);
        myfile << "##fileformat=VCFv4.2\n";
        myfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        myfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
        {myfile << "\t" << *samp;}
        myfile << "\n";

        int snp_order = 0;

        for(vector<string>::iterator it = snp_list.begin();it!=snp_list.end();++it)
        {
            vector <string> data;
            boost::split(data,snp[*it],boost::is_any_of(";"));
            myfile << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t";
            myfile <<".\t.\t.\tGT";
            
            string position = *it;
            
            for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
            {
                
                pair<string, string> key;
                key = make_pair(*samp,*it);
                string snp_data = db[key];
                replace( snp_data.begin(), snp_data.end(), '|', '/');
                vector <string> alleles;
                boost::split(alleles,snp_data,boost::is_any_of("/"));
                
                pair<string, int> key2;
                key2 = make_pair(*samp,snp_order);
                string known_value = knowndb[key2];
                
                if (known_value == "*") {myfile << "\t" << alleles[0] << "/" << alleles[1]; continue;}
                
                
                string knowndata = knownblocks[*samp].substr(1);
                vector<string> known;
                boost::split(known,knowndata,boost::is_any_of(";"));
                vector <string> coord;
                boost::split(coord,known[randblock[*samp]],boost::is_any_of(","));
                
                int start = stoi(coord[0]);
                int end = stoi(coord[0]) + stoi(coord[1]);
                
                
                if ((snp_order >= start) && (snp_order < end))
                {
                    if (known_value == "0") {myfile << "\t" << alleles[0] << "|" << alleles[1]; continue;}
                    if (known_value == "1") {myfile << "\t" << alleles[1] << "|" << alleles[0]; continue;}
                }
                
                myfile << "\t" + alleles[0] + "/" + alleles[1]; continue;
                
            }
            
            myfile << endl;
            snp_order++;
        }
        myfile.close();
    }
    
}











void ReadBeagleResults (string tag, int replicates, int iteration)
{

    map <pair<string,string>,int> final_haplotypes;
 

    for (int a = 1; a <= replicates; a++)
    {
        map <string, string> haplo_sample_h1;
        map <string, string> haplo_sample_h2;
        
        if (v_quiet == 0) {cout << "      > replicate number " + to_string(a) + "\r" << flush;}

        vector <string> samples;
        string file = v_output + "/iteration_" + to_string(iteration) + "/" + tag + "/replicate_" + to_string(a) + ".phased.vcf";
        ifstream inp (file);
        string line;
        while ( getline (inp,line))
        {
            if (line.substr(0,2) == "##") {continue;}
            if (line.substr(0,2) == "#C")
            {
                boost::split(samples,line,boost::is_any_of("\t"));
                continue;
            }
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            for (int b = 9; b < data.size();b++)
            {
                vector <string> alleles;
                boost::split(alleles,data[b],boost::is_any_of("|"));
                haplo_sample_h1[samples[b]] = haplo_sample_h1[samples[b]] + " " + alleles[0];
                haplo_sample_h2[samples[b]] = haplo_sample_h2[samples[b]] + " " + alleles[1];
            }
        }
        inp.close();
        
        for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
        {
            string finalpair = haplo_sample_h1[*samp].substr(1) + "\t" + haplo_sample_h2[*samp].substr(1);
            string finalpair_rev = haplo_sample_h2[*samp].substr(1) + "\t" + haplo_sample_h1[*samp].substr(1);
 
            if (used_pairs.count(finalpair) > 0)
            {
                pair<string, string> key;
                key = make_pair(*samp,finalpair);
                final_haplotypes[key]++;
                continue;
            }
            if (used_pairs.count(finalpair_rev) > 0)
            {
                pair<string, string> key;
                key = make_pair(*samp,finalpair_rev);
                final_haplotypes[key]++;
                continue;
            }
            pair<string, string> key;
            key = make_pair(*samp,finalpair);
            final_haplotypes[key]++;
            used_pairs[finalpair] = 1;
        }
    }
    
    if (tag == "unknown")
    {
        final_haplotypes_vcf.clear();
        final_haplotypes_vcf = final_haplotypes;
    }
}





void CreateKnownSecondRound (string file1, string fileout)
{
//    my @snps = keys %snp;
    
    
    ifstream inp (file1);
    ofstream myfile;
    string outfile = fileout;
    myfile.open (outfile);
    string line;
    while ( getline (inp,line))
    {
        vector <string> data;
        boost::split(data,line,boost::is_any_of(":"));
        string known = "";
        if (find(list_next_round.begin(), list_next_round.end(), data[0]) != list_next_round.end())
        {
            known = string(snp_list.size(),'0');
            myfile << data[0] << ":" << known << endl;
            continue;
        }
       myfile << data[0] << ":" << data[1] << endl;
    }
    inp.close();
    myfile.close();
}




*/
