#include "Analysis/utils/Config.hh"

#include <cstdlib>

std::string Config::workPath = std::string(getenv("ANANAS"))+"/workdir/";
std::string Config::dataPath = Config::workPath+"data/";
std::string Config::rootPath = Config::workPath+"root/";
std::string Config::histPath = Config::workPath+"hist/";
std::string Config::confPath = Config::workPath+"config/";
std::string Config::version  = "";
