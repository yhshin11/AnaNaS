#include "Analysis/selectors/LeptonSelector.hh"

#include <cassert>

#include "Analysis/selectors/ConeSelector.hh"

using namespace std;

// Pt selector
ClassImp( LeptonSelector )

string LeptonSelector::selection[LeptonSelector::kNSel] = 
{ "VeryLoose", "Loose", "Medium", "Tight",
  "VeryLooseID", "LooseID", "MediumID", "TightID",
  "VeryLooseIso", "LooseIso", "MediumIso", "TightIso"};

LeptonSelector::LeptonSelector( unsigned settings, unsigned sel ) : _settings(settings), _sel(sel)
{  

}

bool 
LeptonSelector::accept( const Candidate& cand ) const
{
  CandInfo* info = cand.info();

  bool id;
  if( ! info->getBool( selection[_sel], id ) ) 
    {
      setFlags( cand );
      id = info->getBool( selection[_sel] );
    }

  return id;
}
