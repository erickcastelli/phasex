#include <map>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <mutex>
#include <chrono>


#include "phase-ps.hpp"
#include "vcf.hpp"
#include "external.h"
#include "functions.hpp"
#include "ThreadPool.h"

using namespace std;

mutex mtx;

void main_phase_ps()
{
    
    int check = 1;
    if (! fileExists(v_vcf)) {warnings.push_back("VCF file not informed or not detected"); check=0;}
    if (! fileExists(v_shapeit)) {warnings.push_back("Shapeit4 file not informed or not detected"); check=0;}

    string out = GetStdoutFromCommand(v_shapeit);
    int v_shapeit_version = 0;
    if (out.find("4.1.1") != std::string::npos) {v_shapeit_version = 1;}
    if (out.find("4.1.2") != std::string::npos) {v_shapeit_version = 1;}
    if (out.find("4.1.3") != std::string::npos) {v_shapeit_version = 1;}
    if ((v_shapeit_version == 0) && (out != ""))
    {
        warnings.push_back("This Shapeit version is not compatible with " + Program_name);
        warnings.push_back("Shapeit compatible versions: 4.1.1, 4.1.2, 4.1.3");
        check=0;
    }
    
    
    if (! fileExists(v_beagle)) {warnings.push_back("Beagle JAR file not informed or not detected"); check=0;}
    out = GetStdoutFromCommand("java -jar " + v_beagle);
    int v_beagle_version = 0;
    if (out.find("version 4.1") != std::string::npos) {v_beagle_version = 1;}
    if (out.find("version 5.1") != std::string::npos) {v_beagle_version = 1;}

    if ((v_beagle_version == 0) && (v_beagle != ""))
        {
            warnings.push_back("This Beagle version is not compatible with " + Program_name);
            warnings.push_back("Beagle compatible versions: 4.1, 5.1");
            check=0;
        }


    
    
    
    if (check == 0) {
        screen_message (screen_size, 0, "", 1, 0);
        screen_message (screen_size, 0, "", 1, 0);
        screen_message (screen_size, 0, Program_name + "::PHASE-PS", 1, 0);
        screen_message (screen_size, 0, "", 1, 0);
        screen_message (screen_size, 2, "* Author  : " + Program_author, 1, 0);
        screen_message (screen_size, 2, "* Contact : " + Program_contact, 1, 0);
        screen_message (screen_size, 2, "* Version : " + Program_version, 1, 0);
        screen_message (screen_size, 2, "* Website : " + v_website, 1, 0);
        screen_message (screen_size, 0, "", 1, 0);
        screen_message (screen_size, 2, "Options", 1, 0);
        screen_message (screen_size, 4, "vcf             The input VCF file", 1, 0);
        screen_message (screen_size, 4, "iterations      Maximum number of iterations      (Default: " + to_string(v_iteration) + ")", 1, 0);
        screen_message (screen_size, 4, "replicates      Number of phasing replicates      (Default: " + to_string(v_replicates) + ")", 1, 0);
        screen_message (screen_size, 4, "output          Output folder                     (Default: same as VCF)", 1, 0);
        screen_message (screen_size, 4, "threads         Number of threads                 (Default: " + v_threads + ")", 1, 0);
        screen_message (screen_size, 4, "threshold       Threshold for fixing haplotypes   (Default: " + to_string(v_threshold) + ")", 1, 0);
        screen_message (screen_size, 4, "select          Threshold for selecting haplotype (Default: " + to_string(v_select) + ")", 1, 0);
        screen_message (screen_size, 4, "minfix          Minimum additional fixed samples  (Default: " + to_string(min_fix) + ")" , 1, 0);
        screen_message (screen_size, 4, "scheme          Shapeit scheme                    (Default: see manual)", 1, 0);
        screen_message (screen_size, 4, "shapeit_others  Other parameters for shapeit", 1, 0);
        screen_message (screen_size, 4, "map             Genetic map for human taken from HapMap", 1, 0);
        screen_message (screen_size, 4, "beagle          Beagle JAR file", 1, 0);
        screen_message (screen_size, 4, "shapeit         Shapeit4 bin", 1, 0);
        screen_message (screen_size, 0, "", 1, 0);
        screen_message (screen_size, 4, "--quiet         quiet mode", 1, 0);
        screen_message (screen_size, 4, "--biallelic     phase only biallelic variants", 1, 0);
        screen_message (screen_size, 0, "", 1, 0);
        screen_message (screen_size, 2, "usage: phasex phase-ps vcf=file.vcf", 1, 0);
        screen_message (screen_size, 0, "", 1, 0);
        return;
    }
    
    
    screen_message (screen_size, 0, "", 1, 0);
    screen_message (screen_size, 0, Program_name + "::PHASE-PS", 1, v_quiet);
    
   
    screen_message (screen_size, 2, "Checking VCF structure ...", 2, v_quiet);
    string res = CheckVCF (v_vcf);
    screen_message (screen_size, 2, "Checking VCF structure ... done", 1, v_quiet);

    if (res == "") {return;}
    screen_message (screen_size, 2, res, 1, v_quiet);
    
    
    // Output struct and input structure
    screen_message (screen_size, 2, "Preparing output structure and inputs ...", 2, v_quiet);
    if (v_output == "") {v_output = findfilepath(v_vcf) + "/phasex/";}

    string v_command = "mkdir " + v_output;
    string v_system_out = GetStdoutFromCommand(v_command);

    v_command = "rm -rf " + v_output + "/*";
    v_system_out = GetStdoutFromCommand(v_command);

    v_command = "mkdir " + v_output + "/shapeit/";
    v_system_out = GetStdoutFromCommand(v_command);

    v_command = "mkdir " + v_output + "/source/";
    v_system_out = GetStdoutFromCommand(v_command);

    v_command = "rm " + v_output + "/source/*.vcf";
    v_system_out = GetStdoutFromCommand(v_command);

    v_command = "rm " + v_output + "/source/*.gz";
    v_system_out = GetStdoutFromCommand(v_command);
    
    v_command = "rm " + v_output + "/source/*.tbi";
    v_system_out = GetStdoutFromCommand(v_command);

    
    if (ends_with(v_vcf,"vcf"))
    {
        v_command = "cp " + v_vcf + " " + v_output + "/source/original.vcf";
        v_system_out = GetStdoutFromCommand(v_command);
    }

    if (ends_with(v_vcf,"gz"))
    {
        v_command = "cp " + v_vcf + " " + v_output + "/source/original.vcf.gz";
        v_system_out = GetStdoutFromCommand(v_command);
        v_command = "gunzip " + v_output + "/source/original.vcf.gz";
        v_system_out = GetStdoutFromCommand(v_command);
    }

    v_vcf = v_output + "/source/original.vcf";
    PrepareVcfPS(v_vcf);
    v_command = "bgzip -c " + v_output + "/source/biallelic.vcf > " + v_output + "/source/biallelic.vcf.gz";
    v_system_out = GetStdoutFromCommand(v_command);
    v_command = "tabix -p vcf " + v_output + "/source/biallelic.vcf.gz";
    v_system_out = GetStdoutFromCommand(v_command);
    
    screen_message (screen_size, 2, "Preparing output structure and inputs ... done", 1, v_quiet);

    
    
    
    
    
    
    
    
    
    ofstream log;
    string logfile = v_output + "/phasex.log";
    log.open (logfile);
    
    v_message = Program_name + "::PHASE-PS, version " + Program_version;
    log << v_message << endl;
    log << "output folder: " << v_output << endl;
    
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    
    strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);
    std::string str(buffer);
    
    log << "Date and time: " << str << endl;
    log << "VCF file: " << v_vcf << endl;;
    log << "Iterations: " << to_string(v_iteration) << endl;
    log << "Replicates: " << to_string(v_replicates) << endl;
    log << "Beagle JAR file: " << v_beagle << endl;
    log << "Shapeit file: " << v_shapeit << endl;
    log << "Threshold for fixing an haplotype: " << to_string(v_threshold) << endl;
    log << "Threshold for selecting an haplotype pair: " << to_string(v_select) << endl;
    log << "Minfix: " << to_string(min_fix) << endl;
    log << "Output folder: " << v_output << endl;
    log << "Region: " << v_region << endl;
    log << "Map: " << v_map << endl;
    log << "Threads: " << v_threads << endl;
    if (onlybiallelic == 1) {log << "Type: biallelic" << endl;}
    if ((onlybiallelic == 0) && (multiallelic == 1)) {log << "Type: multi-allelic" << endl;}
    
    log << endl;

    
    
    
    // ********************* SHAPEIT *************************
    
    
    int iterations_done = 0;
    for (int iter = 1; iter <= v_iteration; iter++)
    {
        
        iterations_done++;
        v_command = "mkdir " + v_output + "/shapeit/iteration_" + to_string(iter);
        v_system_out = GetStdoutFromCommand(v_command);
        
        v_command = "rm " + v_output + "/shapeit/iteration_" + to_string(iter) + "/*.*";
        v_system_out = GetStdoutFromCommand(v_command);
        

        // starting replications for shapeit

        string input = "";
        if (iter == 1) {input= v_output + "/source/biallelic.vcf.gz";}
        if (iter > 1) {input= v_output + "/shapeit/iteration_" + to_string((iter-1)) +  "/phased_final.vcf.gz";}
        
        if (! fileExists(input)) {warnings.push_back("Something went wrong ... no input file for this iteration."); return;}

        
        string v_message = "Shapeit: iteration " + to_string(iter) + ", starting replicates  ...";
        screen_message (screen_size, 2, v_message, 2, v_quiet);
        ThreadPool pool(stoi(v_threads));
        std::vector< std::future<int> > results;
        ready = 0;
        for (int replica = 1; replica <= v_replicates; replica++)
        {
            results.emplace_back(
                pool.enqueue([replica,iter,input] {
                    int seed = int(rand());
                    string v_outfile = v_output + "/shapeit/iteration_" + to_string(iter) + "/phased_" + to_string(replica) + ".vcf";
                    string v_logfile = v_output + "/shapeit/iteration_" + to_string(iter) + "/phased_" + to_string(replica) + ".log";
                    string v_command = v_shapeit + " --input " + input + " --use-PS " + shapeit_useps + " --output " + v_outfile + " --region " + v_region + " --mcmc-iterations " + shapeit_scheme + " --seed " + to_string(seed) + " --thread 1 --log " + v_logfile + " --pbwt-depth " + shapeit_pbwt_depth + " " + shapeit_other_parameters + " --sequencing";
                    if (v_map != "") {v_command = v_command + " --map " + v_map;}
                    string v_system_out = GetStdoutFromCommand(v_command);
                    mtx.lock();
                    ready++;
                    string v_message = "Shapeit: iteration " + to_string(iter) + ", replicate " + to_string(ready) + " out of " + to_string(v_replicates) + " ...";
                    screen_message (screen_size, 2, v_message, 2, v_quiet);
                    mtx.unlock();
                    return replica*replica;
                })
            );

        } //end replicate shapeit

        extern string shapeit_scheme;
        extern string shapeit_useps;
        extern string shapeit_pbwt_depth;
        
        for(auto && result: results){result.get();} // waiting for all threads
        
        
        
        
        
        v_message = "Shapeit: iteration " + to_string(iter) + ": reading and comparing data ...";
        screen_message (screen_size, 2, v_message, 2, v_quiet);
        fixed_haplotypes.clear(); //new
        LoadingGenotypes(v_output + "/source/biallelic.vcf");
        ReadResults ("shapeit",v_replicates, iter);
        
        int sample_size = sample_list.size();
        list_next_round.clear();
        
        
        ofstream myfile;
        string outfile = v_output + "/shapeit/iteration_" + to_string(iter) + "/haplotypes.txt";
        myfile.open (outfile);
        
        myfile << "Sample\th1\th2\tfreq\tinfo" << endl << endl;
        
        for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
        {
            map <string,int> haplos_for_this_sample;
            for(std::map<string,int>::iterator iter = used_pairs.begin(); iter != used_pairs.end(); ++iter)
            {
                string haplo =  iter->first;
                pair<string, string> key;
                key = make_pair(*samp,haplo);
                if (final_haplotypes_vcf.count(key) > 0) {haplos_for_this_sample[haplo] = 1;}
            }
            
            for(std::map<string,int>::iterator iter = haplos_for_this_sample.begin(); iter != haplos_for_this_sample.end(); ++iter)
            {
                string haplo =  iter->first;
                float vcf_p = 0;
                pair<string, string> key;
                key = make_pair(*samp,haplo);
                
                myfile << *samp << "\t" << haplo;
                
                if(final_haplotypes_vcf.count(key) > 0)
                {
                    vcf_p = (float)final_haplotypes_vcf[key] / (float)v_replicates;
                    myfile << "\t" + to_string(vcf_p);
                }
                else {myfile << "\t-";}
                int pass = 1;
                if (vcf_p < v_threshold) {pass = 0;}
                
                
                if (pass == 1)
                {
                    myfile << "\tdef";
                    list_next_round.push_back(*samp);
                    fixed_haplotypes[*samp] = haplo;
                }
                if (pass == 0) {myfile << "\t-";}
                myfile << endl;
            }
            
            myfile << endl;
        }
        
        myfile.close();
        v_message = "Shapeit: iteration " + to_string(iter) + ": samples with defined haplotypes: " + to_string(list_next_round.size()) + " (" + to_string((((float)list_next_round.size() / (float)sample_size)*100)) + "%)";
        screen_message (screen_size, 2, v_message, 1, v_quiet);
        
        log << "Shapeit: Iteration " + to_string(iter) + ": samples with defined haplotypes: " << to_string(list_next_round.size()) + " out of " + to_string(sample_size) + " (" + to_string((((float)list_next_round.size() / (float)sample_size)*100)) + "%)" << endl;
        
        if ((iter > 1) && (last_pass >= list_next_round.size())) {break;}
        if ((iter > 1) && (list_next_round.size() < last_pass + min_fix)) {break;}

        last_pass = list_next_round.size();
        
        if (iter < v_iteration) {
            v_message = "Shapeit: starting iteration " + to_string(iter+1) + ": preparing files ...";
            screen_message (screen_size, 2, v_message, 2, v_quiet);
        }
        LoadingGenotypesFixed();
        PrintVCFPSnew ("shapeit", iter);
        
        
    } // end iter

    
    
    
    
    v_message = "Shapeit: printing final haplotypes ... ";
    screen_message (screen_size, 2, v_message, 2, v_quiet);

    map <pair<string,int>,string> finaldata;
    map <string,string> haplos_per_samples;
    unordered_map <string,string> sample_haplotypes;
    
    for (int a = 1; a <= iterations_done; a++)
    {
        
        ifstream inp (v_output + "/shapeit/iteration_" + to_string(a) + "/haplotypes.txt");
        string line;
        getline (inp,line);
        while ( getline (inp,line))
        {
            if (line == "") {continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            pair<string, int> key;
            key = make_pair(data[0] + "\t" + data[1] + "\t" + data[2], a);
            finaldata[key] = data[3] + "\t" + data[4];
            pair<string, string> key2;
            haplos_per_samples[data[0]] = haplos_per_samples[data[0]] + data[1] + "\t" + data[2] + "\n";
        }
        inp.close();
    }
    
    
    ofstream outresults;
    string outfile = v_output + "/shapeit/results.txt";
    outresults.open (outfile);
    
    outresults << "Sample\th1\th2";
    for (int a = 1; a <= iterations_done; a++)
    {
        outresults << "\tFreq(" + to_string(a) + ")";
        outresults << "\tInfo(" + to_string(a) + ")";
    }
    outresults << "\tStatus" << endl << endl;

    
    passed_samples.clear();
    
    for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
    {
        vector <string> haplos;
        boost::split(haplos,haplos_per_samples[*samp],boost::is_any_of("\n"));
        unordered_map <string,int> list;
        
        for (auto &haplodata : haplos)
        {
            if (haplodata == "") {continue;}
            list[haplodata] = 1;
        }
        
        
        for (auto &haplodata : list)
        {
            
            outresults << *samp << "\t" << haplodata.first;
            
            for (int b = 1; b <= iterations_done; b++)
            {
                pair<string, int> key;
                key = make_pair(*samp + "\t" + haplodata.first, b);
                if (finaldata.find(key) == finaldata.end())
                {
                    finaldata[key] = "-\t-";
                    outresults << "\t-\t-";
                }
                else
                {
                    outresults << "\t" << finaldata[key];
                }
                
                
                if (b == iterations_done)
                {
                    string res = "";
                    vector <string> p;
                    boost::split(p,finaldata[key],boost::is_any_of("\t"));
                    if (p[0] == "-") {res = "-";}
                    if (p[0] == "") {res = "-";}
                    if (p[0] != "-") {if(stof(p[0]) >= v_select) {res = "pass";passed_samples.push_back(*samp);}}
                    outresults << "\t" << res;
                }
            }
            outresults << endl;
        }
        outresults << endl;
    }
    outresults.close();


    v_message = "Shapeit: printing final haplotypes ... done";
    screen_message (screen_size, 2, v_message, 1, v_quiet);
    
    int sample_size = sample_list.size();
    v_message = "Shapeit: samples with haplotypes above threshold: " + to_string(passed_samples.size()) + " (" + to_string((((float)passed_samples.size() / (float)sample_size)*100)) + "%)";
    screen_message (screen_size, 2, v_message, 1, v_quiet);
    log << v_message << endl;

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   
    
    
    
    
    
    // Beagle for multi-allelic variants
    if ((onlybiallelic == 0) && (multiallelic == 1))
    {
        v_message = "Beagle:  starting phasing multi-allelic variants: preparing files ... ";
        screen_message (screen_size, 2, v_message, 2, v_quiet);
        v_command = "mkdir " + v_output + "/beagle/";
        v_system_out = GetStdoutFromCommand(v_command);
        v_command = "mkdir " + v_output + "/beagle/source/";
        v_system_out = GetStdoutFromCommand(v_command);
        v_command = "mkdir " + v_output + "/beagle/iteration_1/";
        v_system_out = GetStdoutFromCommand(v_command);
        
        LoadingGenotypesFixed();
        LoadingGenotypesMulti(v_output + "/source/multiallelic.vcf");
        PrintVCFBeagle ();
        if (! fileExists(v_output + "/beagle/source/source.vcf")) {warnings.push_back("Something went wrong ... no input file for this iteration."); return;}
        LoadingGenotypes(v_output + "/beagle/source/source.vcf");
        fixed_haplotypes.clear(); //new
        list_next_round.clear(); //new
        passed_samples.clear(); //new
        finaldata.clear(); //new
        last_pass = 0; //new
        list_next_round.clear(); //new
        
        v_message = "Beagle:  starting phasing multi-allelic variants: preparing files ... done";
        screen_message (screen_size, 2, v_message, 1, v_quiet);

        v_message = "Beagle:  starting iterations ... ";
        screen_message (screen_size, 2, v_message, 2, v_quiet);

 
        
        
        int iterations_done = 0;
        for (int iter = 1; iter <= v_iteration; iter++)
        {
            
            iterations_done++;
            v_command = "mkdir " + v_output + "/beagle/iteration_" + to_string(iter);
            v_system_out = GetStdoutFromCommand(v_command);
            string input;
            if (iter == 1) {input = v_output + "/beagle/source/source.vcf";}
            if (iter > 1) {input= v_output + "/beagle/iteration_" + to_string((iter-1)) +  "/phased_final.vcf";}

            if (! fileExists(input)) {warnings.push_back("Something went wrong ... no input file for this iteration."); return;}
            
            
            string v_message = "Beagle:  iteration " + to_string(iter) + ", starting replicates  ...";
            screen_message (screen_size, 2, v_message, 2, v_quiet);
            ThreadPool pool(1);
            std::vector< std::future<int> > results;
            ready = 0;
            
            for (int replica = 1; replica <= v_replicates; replica++)
            {
                results.emplace_back(
                    pool.enqueue([replica,iter,input] {
                    int seed = (replica * 1000) + replica - 5;
                    string v_outfile = v_output + "/beagle/iteration_" + to_string(iter) + "/phased_" + to_string(replica);
                    string v_logfile = v_output + "/beagle/iteration_" + to_string(iter) + "/phased_" + to_string(replica) + ".log";
                    
                    //retorno:
                    string v_command = "java -Xmx8g -jar " + v_beagle + " seed=" + to_string(seed) + " gt=" + input + " out=" + v_outfile + " nthreads=" + v_threads + " > " + v_logfile;
                    string v_system_out = GetStdoutFromCommand(v_command);
                   
                    //if (! fileExists(v_outfile + ".vcf.gz")) {goto retorno;}
                    v_command = "gunzip -f " + v_outfile + ".vcf.gz";
                    v_system_out = GetStdoutFromCommand(v_command);
                    mtx.lock();
                    ready++;
                    string v_message = "Beagle:  iteration " + to_string(iter) + ", replicate " + to_string(ready) + " out of " + to_string(v_replicates) + " ...";
                    screen_message (screen_size, 2, v_message, 2, v_quiet);
                    mtx.unlock();
                    return replica*replica;
                })
                );
                
            } //end replicate
            
            for(auto && result: results){result.get();} // waiting for all threads
            
            
            
            v_message = "Beagle:  iteration " + to_string(iter) + ": reading and comparing data ...";
            screen_message (screen_size, 2, v_message, 2, v_quiet);
            LoadingGenotypes(v_output + "/source/original-ps.vcf");
            ReadResults ("beagle",v_replicates, iter);

            int sample_size = sample_list.size();
            list_next_round.clear();
            
            ofstream myfile;
            string outfile = v_output + "/beagle/iteration_" + to_string(iter) + "/haplotypes.txt";
            myfile.open (outfile);
            
            myfile << "Sample\th1\th2\tfreq\tinfo" << endl << endl;
            
            for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
            {
                map <string,int> haplos_for_this_sample;
                for(std::map<string,int>::iterator iter = used_pairs.begin(); iter != used_pairs.end(); ++iter)
                {
                    string haplo =  iter->first;
                    pair<string, string> key;
                    key = make_pair(*samp,haplo);
                    if (final_haplotypes_vcf.count(key) > 0) {haplos_for_this_sample[haplo] = 1;}
                }
                
                for(std::map<string,int>::iterator iter = haplos_for_this_sample.begin(); iter != haplos_for_this_sample.end(); ++iter)
                {
                    string haplo =  iter->first;
                    float vcf_p = 0;
                    pair<string, string> key;
                    key = make_pair(*samp,haplo);
                    
                    myfile << *samp << "\t" << haplo;
                    
                    if(final_haplotypes_vcf.count(key) > 0)
                    {
                        vcf_p = (float)final_haplotypes_vcf[key] / (float)v_replicates;
                        myfile << "\t" + to_string(vcf_p);
                    }
                    else {myfile << "\t-";}
                    int pass = 1;
                    if (vcf_p < v_threshold) {pass = 0;}
                    
                    
                    if (pass == 1)
                    {
                        myfile << "\tdef";
                        list_next_round.push_back(*samp);
                        fixed_haplotypes[*samp] = haplo;
                    }
                    if (pass == 0) {myfile << "\t-";}
                    myfile << endl;
                }
                
                myfile << endl;
            }
            
            myfile.close();
            v_message = "Beagle:  iteration " + to_string(iter) + ": " + to_string(list_next_round.size()) + " samples above threshold  (" + to_string((((float)list_next_round.size() / (float)sample_size)*100)) + "%)";
            screen_message (screen_size, 2, v_message, 1, v_quiet);
            
            log << "Beagle: Iteration " + to_string(iter) + ": samples with defined haplotypes: " << to_string(list_next_round.size()) + " out of " + to_string(sample_size) + " (" + to_string((((float)list_next_round.size() / (float)sample_size)*100)) + "%)" << endl;
            
            
            if ((iter > 1) && (last_pass >= list_next_round.size())) {break;}
            if ((iter > 1) && (list_next_round.size() < last_pass + min_fix)) {break;}

            last_pass = list_next_round.size();
            
            v_message = "Beagle:  starting iteration " + to_string(iter+1) + ": preparing files ...";
            screen_message (screen_size, 2, v_message, 2, v_quiet);

            print_debug (debug, "Beagle, Starting LoadingGenotypeFixed iteration " + to_string(iter));
            LoadingGenotypesFixed();
            print_debug (debug, "Beagle, Ending LoadingGenotypeFixed iteration " + to_string(iter));
 
            PrintVCFPS ("beagle", iter);
            
        } //end iter
        
        
        
        v_message = "Beagle:  printing final haplotypes ... ";
        screen_message (screen_size, 2, v_message, 2, v_quiet);
        
        map <pair<string,string>,string> finaldata;
        unordered_map <string,string> sample_haplotypes;
        
        for (int a = 1; a <= iterations_done; a++)
        {
            
            ifstream inp (v_output + "/beagle/iteration_" + to_string(a) + "/haplotypes.txt");
            string line;
            while ( getline (inp,line))
            {
                if (line == "") {continue;}
                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));
                
                pair<string, string> key;
                key = make_pair(data[0],data[1] + "\t" + data[2]);
                finaldata[key] = finaldata[key] + to_string(a);
                for (int b = 3; b < data.size(); b++) {finaldata[key] = finaldata[key] + "," + data[b];}
                finaldata[key] = finaldata[key] + "\n";
            }
            inp.close();
        }
        
        
        ofstream outresults;
        string outfile = v_output + "/beagle/results.txt";
        outresults.open (outfile);
        
        outresults << "Sample\th1\th2";
        for (int a = 1; a <= iterations_done; a++)
        {
            outresults << "\tFreq(" + to_string(a) + ")";
            outresults << "\tInfo(" + to_string(a) + ")";
        }
        outresults << "\tStatus" << endl << endl;
        
        
        
        vector <string> passed_samples;
        
        for(vector<string>::iterator samp = sample_list.begin();samp!=sample_list.end();++samp)
        {
            for (auto &p : finaldata)
            {
                if (p.first.first != *samp){continue;}
                outresults << *samp << "\t";
                outresults << p.first.second;
                
                pair<string, string> key;
                key = make_pair(*samp,p.first.second);
                string data = finaldata[key];
                
                vector <string> data_iter;
                boost::split(data_iter,data,boost::is_any_of("\n"));
                
                map <int,string> sub;
                int number_of_fields = 0;
                for(vector<string>::iterator it = data_iter.begin();it!=data_iter.end();++it)
                {
                    if (*it == "") {continue;}
                    vector <string> values;
                    boost::split(values,*it,boost::is_any_of(","));
                    number_of_fields = values.size();
                    for (int c = 1; c < number_of_fields; c++)
                    {
                        sub[stoi(values[0])] = sub[stoi(values[0])] + "\t" + values[c];
                    }
                    
                }
                
                
                for (int a = 1; a <= iterations_done; a++)
                {
                    if (sub.count(a) > 0 )
                    {
                        outresults << sub[a];
                    }
                    
                    if (sub.count(a) == 0 )
                    {
                        for (int d = 1; d < number_of_fields; d++)
                        {
                            outresults << "\t-";
                        }
                    }
                }
                
                data = sub[iterations_done];
                vector <string> values;
                boost::split(values,data,boost::is_any_of("\t"));
                
                //cout << data << endl;
                
                //cout << ">" << data << "<" << to_string(values.size()) << endl;
                
                int pass = 1;
                if (values.size() < 2) {pass=0;}
                if (values.size() > 2) {
                    for (int c = 1; c < (values.size() - 1); c++)
                    {
                        if (values[c] == "-") {pass = 0;continue;}
                        if (values[c] == "") {pass = 0;continue;}
                        if ( stof(values[c]) < v_select ){pass = 0;}
                    }
                    if (pass == 1) {outresults << "\tpass";passed_samples.push_back(*samp);} else {outresults << "\t-";}
                }
                
                outresults << endl;
            }
            
            outresults << endl;
        }
        
        outresults.close();
        
        
        v_message = "Beagle:  printing final haplotypes ... done";
        screen_message (screen_size, 2, v_message, 1, v_quiet);
        
        int sample_size = sample_list.size();
        v_message = "Beagle:  samples with haplotypes above threshold: " + to_string(passed_samples.size()) + " (" + to_string((((float)passed_samples.size() / (float)sample_size)*100)) + "%)";
        screen_message (screen_size, 2, v_message, 1, v_quiet);
        log << v_message << endl;
        
        
    }
    

    v_input = v_output;
    main_vcf();
    
    screen_message (screen_size, 2, "", 1, v_quiet);
    log.close();

    return;

  
    
}
