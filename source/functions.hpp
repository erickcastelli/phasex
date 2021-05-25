#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sys/stat.h>


#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>

using namespace std;
void print_debug (int debug, string msg);
void print_warnings ();
void screen_message (int size, int left, string message, int enter, int quiet);
string findfilepath (string v_file);
bool fileExists(const std::string& filename);
string GetStdoutFromCommand(string cmd);
void printmessage (string v_text, int v_quiet, int v_level, int v_enter);
//void CreateSubFolderStructure (string output, int iteration);
bool ends_with(const std::string &filename, const std::string &ext);
void ReadResults (string tag, int replicates, int iteration);
void LoadingGenotypes (string file);
void LoadingGenotypesMulti (string file);
void PrepareVcfPS (string file);
void PrintVCFPS(string tag, int iter);
void PrintVCFPSnew(string tag, int iter);
void PrintVCFBeagle();
string CheckVCF (string file);
void LoadingGenotypesFixed ();


//string update(string url, string app, int ignore_update);
#endif /* functions_hpp */
