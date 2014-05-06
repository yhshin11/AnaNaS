#include "Analysis/utils/CfgParser.hh"

#include <cassert>
#include <iostream>

const std::string CfgParser::commentSeparator  = "#";
const std::string CfgParser::keySeparatorLeft  = "[";
const std::string CfgParser::keySeparatorRight = "]";
const std::string CfgParser::defaultAttribute  = "*";

CfgParser::CfgParser( const char * fileName )
{
    open( fileName );
    stripped_ = false;
}

void CfgParser::open( const char * fileName )
{
    fileName_ = fileName;
    ifile_.open( fileName ); 
    if ( ifile_.is_open() ) {
        parse();
    } else {
        std::cerr << "Error opening file " << fileName_ << std::endl;
        exit(1);
    }
}

void CfgParser::strip( std::string & s )
{
    if ( stripped_ ) return;
    size_t found;
    // comments
    found = s.find( commentSeparator );
    if ( found != std::string::npos  ) {
        s.replace( found, s.length(), "" );
    }
    // blanks
    while ( 1 ) {
        found = s.find( " " );
        if ( found != std::string::npos  ) {
            s.replace( found, 1, "" );
        } else {
            break;
        }
    }
    // tabs
    while ( 1 ) {
        found = s.find( "\t" );
        if ( found != std::string::npos  ) {
            s.replace( found, 1, "" );
        } else {
            break;
        }
    }
    stripped_ = true;
}

bool CfgParser::stripDefault( std::string & s )
{
    size_t found = 0;
    while ( 1 ) {
        found = s.find( "*", found );
        if ( found == std::string::npos ) return false;
        if ( found != s.length()-1 ) {
            std::cerr << "Warning: value string " << s
                << " contains default attribute not at the end.\n";
        } else {
            s.replace( found, 1, "" );
            return true;
        }
    }
    return false;
}

void CfgParser::key( std::string & s, std::string & k )
{
    strip( s );
    size_t foundLeft, foundRight;
    foundLeft = s.find( keySeparatorLeft );
    if ( foundLeft != std::string::npos ) {
        foundRight = s.find( keySeparatorRight );
        if ( foundRight != std::string::npos
             && foundLeft < foundRight ) {
            k.assign(s, foundLeft+1, foundRight-foundLeft-1);
        }
    }
}

void CfgParser::value( std::string & s, std::string & v )
{
    strip( s );
    size_t foundLeft, foundRight;
    foundLeft = s.find( keySeparatorLeft );
    foundRight = s.find( keySeparatorRight );
    if ( foundLeft == std::string::npos && foundRight == std::string::npos ) {
        // no forbidden characters in the string
        v = s;
    } else {
        std::cerr << "Field " << s << " does not apperear to be"
            " neither a key nor a value.\n";
    }
}

void CfgParser::parse()
{
    std::string line;
    while ( 1 ) {
        std::getline( ifile_, line );
        stripped_ = false;
        if ( ifile_.eof() ) break;
        strip( line );
        std::string k;
        key( line, k );
        // check if it is a key
        if ( k.length() ) {
            currentKey_ = k;
            cfg_[ currentKey_ ];
        } else {
            // not a key: may be a value
            std::string v;
            value( line, v );
            if ( v.length() ) {
                // not a comment (already stripped)
                if ( currentKey_.length() ) {
                    // if it is default, get it
                    bool isDefault = stripDefault( v );
                    cfg_[ currentKey_ ].push_back( v );
                    if ( isDefault ) {
                        Default::const_iterator it = default_.find( currentKey_ );
                        if ( it == default_.end() ) {
                            default_[ currentKey_ ] = cfg_[ currentKey_ ].size()-1;
                        } else {
                            std::cerr << "Already made a default for the key " << currentKey_
                                << " value " << v
                                << ". Please check your config file syntax.\n";
                        }
                    }
                } else {
                    // not a value
                    std::cerr << "No key for value " << v << ". "
                        "Please check the syntax of your file\n";
                }
            }
        }
    }
}

const CfgParser::Values & CfgParser::values( Key & k ) const
{
    Config::const_iterator it = cfg_.find( k );
    assert( it != cfg_.end() );
    return it->second;
}

void CfgParser::values( Key & k, Values & v ) const
{
    Config::const_iterator it = cfg_.find( k );
    if( it == cfg_.end() ){
      std::cerr << __FILE__ << ":" << __LINE__ << ": "
                << "Paramater " << k
                << " not found in configuration file "
                << fileName_ << std::endl;
      exit(1);
    }
    v = it->second;
}

void CfgParser::values( const char * k, Values & v ) const
{
    std::string kk(k);
    values( kk, v );
}

size_t CfgParser::defaultValue( Key & k ) const
{
    Default::const_iterator it = default_.find( k );
    if ( it == default_.end() ) {
        std::cerr << "No default value found for key " << k << ".\n";
        if ( values( k ).size() ) {
            std::cerr << "Using first value as default.\n";
            return 0;
        }
    } else {
        return it->second;
    }
    // it never gets here
    return 0;
}

size_t CfgParser::defaultValue( const char * k ) const
{
    std::string kk(k);
    return defaultValue( kk );
}

//GHM
bool CfgParser::hasKey( const char* k ) const
{
  Config::const_iterator kit = cfg_.find( k );
  if( kit==cfg_.end() ) return false;
  return true;
}
//GHM

void CfgParser::print( std::ostream &o )
{
    for ( Config::const_iterator kit = cfg_.begin(); kit != cfg_.end(); ++kit ) {
        Key k = kit->first;
        o << k << "\n";
        Values v = kit->second;

	//GHM modified by GHM to indicate which is the default
	size_t idef_ = v.size();
	Default::const_iterator dit_ = default_.find( k );
	if( dit_!=default_.end() ) idef_ = dit_->second;
        for( size_t ii=0; ii<v.size(); ii++ ) 
	  {
	    if( ii==idef_ ) o << "==> ";	    
	    else            o << "... ";
	    o << v[ii] << "\n";
	  }
	//GHM
    }
    //GHM
    //GHMfor ( Default::const_iterator kit = default_.begin(); kit != default_.end(); ++kit ) {
    //GHM    o << "__ " << kit->first << " " << kit->second << "\n";
    //GHM    }
}

int main()
{
    CfgParser p("cfg.log");
    p.print( std::cout );
    CfgParser::Key k = "provissima";
    CfgParser::Values v = p.values( k );
    std::cout << "--> " << v.size() << "\n";
    for ( size_t s = 0; s < v.size(); ++s ) {
        std::cout << v[s] << "\n";
    }
    return 0;
}
