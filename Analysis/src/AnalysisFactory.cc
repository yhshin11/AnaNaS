#include "Analysis/src/AnalysisFactory.hh"

using namespace std;

#include "Analysis/src/DiphotonAnalysis.hh"
#include "Analysis/src/ZllAnalysis.hh"
#include "Analysis/src/DibosonAnalysis.hh"
#include "Analysis/src/Commissioning.hh"
#include "Analysis/src/WAnalysis.hh"
#include "Analysis/src/OfficialWAnalysis.hh"
#include "Analysis/src/OfficialZAnalysis.hh"
#include "Analysis/src/DYContamination.hh"
#include "Analysis/src/MCAnalysis.hh"
#include "Analysis/src/METAnalysis.hh"
#include "Analysis/src/JETAnalysis.hh"
#include "Analysis/src/DMAnalysis.hh"

SampleAnalysis*
AnalysisFactory::get( const string& analysis, Sample& sample, EventManager& eventManager )
{
  string analysis_(analysis);

  if( analysis_==string("Sample") )
    {
      return new SampleAnalysis( "Sample", sample, eventManager );
    }
  if( analysis_==string("Zll") )
    {
      //      ZllAnalysis::imode = EventServer::kZee;
      return new ZllAnalysis( sample, eventManager );
    }
  //   if( analysis_==string("Zmm") )
  //     {
  //       //       ZllAnalysis::imode = EventServer::kZmm;
  //       return new ZllAnalysis( sample, eventManager );
  //     }
  if( analysis_==string("Diboson") )
    {
      return new DibosonAnalysis( sample, eventManager );
    }
  if( analysis_==string("MET") )
    {
      return new METAnalysis( sample, eventManager );
    }
  if( analysis_==string("JET") )
    {
      return new JETAnalysis( sample, eventManager );
    }
  if( analysis_==string("DM") )
    {
      return new DMAnalysis( sample, eventManager );
    }


//   if( analysis_==string("Diphoton") )
//     {
//       return new DiphotonAnalysis( sample, eventManager );
//     }
//   if( analysis_==string("Commissioning") )
//     {
//       return new Commissioning( sample, eventManager );
//     }
//   if( analysis_==string("mmWe") )
//     {
//       WAnalysis::imode = EventServer::kWe;
//       return new WAnalysis( sample, eventManager );
//     }
//   if( analysis_==string("mmWm") )
//     {
//       WAnalysis::imode = EventServer::kWm;
//       return new WAnalysis( sample, eventManager );
//     }
//   if( analysis_==string("OWe") )
//     {
//       OfficialWAnalysis::imode = EventServer::kWe;
//       return new OfficialWAnalysis( sample, eventManager );
//     }
//   if( analysis_==string("OWm") )
//     {
//       OfficialWAnalysis::imode = EventServer::kWm;
//       return new OfficialWAnalysis( sample, eventManager );
//     }
//   if( analysis_==string("OZee") )
//     {
//       OfficialZAnalysis::imode = EventServer::kZee;
//       return new OfficialZAnalysis( sample, eventManager );
//     }
//   if( analysis_==string("OZmm") )
//     {
//       OfficialZAnalysis::imode = EventServer::kZmm;
//       return new OfficialZAnalysis( sample, eventManager );
//     }
//   if( analysis_==string("DYContamination") )
//     {
//       return new DYContamination( sample, eventManager );
//     }
//   if( analysis_==string("MC") )
//     {
//       return new MCAnalysis( sample, eventManager );
//     }

  return 0;
}
