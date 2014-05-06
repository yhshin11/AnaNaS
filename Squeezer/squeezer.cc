// Squeezer application, a utility to copy and skim AnaNaS n-tuple files
// Author: Ph. Gras CEA/IRFU Saclay
// Jan. 3, 11

#include "CopyAnaTuple.h"

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <TFile.h>

using namespace std;

void help(){
  cout << "\n"
    "Usage: 1) squeezer [OPTIONS] INPUT_FILE OUTPUT_FILE\n"
    "       2) squeezer {--make-event-list|--make-branch-list} LIST_OUTPUT_FILE INPUT_FILE\n"
    "       3) squeezer --help\n"
    "\n"
    "Decription: produces skim of AnaNaS n-tuple file.\n"
    "\n"
    "Note: the class CopyAnaTuple can be extended to customize event and tree branch selection.\n"
    "\n"
    "In form (1) events of INPUT_FILE are copied into OUTPUT_FILE. Use options to selects the events.\n"
    "In form (2) produces a file with either:\n"
    "  the list of events in the format read by --event-list-from option\n"
    "  or the list of branch in the format read by --exclude-branch-from.\n"
    "Use - in place of "
    "EVENT_LIST_OUTPUT_FILE to display the list on standard ouput.\n"
    "In form (3) displays this help.\n"
    "\n"
    "OPTIONS:\n"
    "-v                           increase verbosity. Add more v for more verbosity, e.g. -vv.\n"
    "--max-events NEVENTS         limit processing to the first NEVENTS events of the INPUT_FILE file.\n"
    "--event-list-from FILE       read list of events to copy from file FILE. Each line \n"
    "must contain a run number followed by an event number. Line starting with a '#' sign\n"
    "are considered as comments and are ignored.\n"
    "--exclude-branches-from FILE excludes branches list in FILE from the copy. Format is one branch "
    "per line. Each line must contain the tree name followed by either the branch name or *.\n"
    "In the latter case the full tree is excluded from copy.\n"
    "Line starting with a '#' sign are considered as comments and are ignored.\n"
       << std::endl;
}

int parse_cmd_line(CopyAnaTuple& cat, int argc, char* argv[]);

int main(int argc, char* argv[]){

  CopyAnaTuple cat;
  
  int i = parse_cmd_line(cat, argc, argv);
  argc -= i;
  argv += i;

  if(argc!=2){
    cerr << "Wrong command usage. See copyAnatuple --help." << endl;
    return 1;
  }

  cat.run(argv[0], argv[1]);

  return 0;
}

int parse_cmd_line(CopyAnaTuple& cat, int argc, char* argv[]){
   typedef enum {no_arg=0, required_arg, optional_arg} has_arg_t;
   enum {make_event_list = 300, make_branch_list, event_list_from,
         exclude_branches_from};
   static struct option options[] = {
     {"verbose", no_arg, NULL, 'v'},
     {"help", no_arg, NULL, 'h'},
     {"max-events", required_arg, NULL, 'n'},
     {"make-event-list", required_arg, NULL, make_event_list},
     {"make-branch-list", required_arg, NULL, make_branch_list},
     {"event-list-from", required_arg, NULL, event_list_from},
     {"exclude-branches-from", required_arg, NULL, exclude_branches_from},
     {0, 0, 0, 0}
   };

   const int noptions = sizeof(options)/sizeof(options[0])-1;
   char short_options[3*noptions+1];

   //for each long option, use as short option equivalent the option.val
   //character
   int pos = 0;
   for(int ioption=0; ioption<noptions; ++ioption){
     short_options[pos++] = options[ioption].val;
     switch(options[ioption].has_arg){
     case required_argument: //a required arg. is indicated with a colon
       short_options[pos++] = ':';
       break;
     case optional_argument://an optional arg is indicated with two colons
       short_options[pos++] = ':';
       short_options[pos++] = ':';
       break;
     }
   }
   short_options[pos] = '\0';

   //   cout << "short option description: " << short_options << endl;

   int c = 0;
   while((c=getopt_long(argc, argv, short_options,
                        options, NULL))!=-1){
     switch (c){
     case 'h'://-h or --help
       help();
       exit(0);
     case 'v'://-v or --verbose
       ++cat.verbose_;
       break;
     case 'n':
       {
         int val = strtol(optarg, 0, 0);
         cat.setMaxEvents(val);
         break;
       }
     case make_event_list:
       if(strcmp(optarg, "-")==0){
         cat.listEvents(cout, optarg);
       } else{
         ofstream o(argv[2]);
         cat.listEvents(o, argv[3]);
       }
       exit(0);
     case make_branch_list:
       if(strcmp(optarg, "-")==0){
         cat.listBranches(cout, argv[optind]);
       } else{
         ofstream o(optarg);
         cat.listBranches(o, argv[optind]);
       }
       exit(0);
     case event_list_from:
       cat.readEventList(optarg);
       break;     
     case exclude_branches_from:
       cat.readExludedBranchList(optarg);
       break;
     }
   }
   if(c=='?') exit(1);
   return optind;
}
