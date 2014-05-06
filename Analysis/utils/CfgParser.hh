#ifndef CfgParser_h
#define CfgParser_h

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>

class CfgParser {
    public:

        typedef std::string Key;
        typedef std::string Val;
        typedef std::vector< Val > Values;
        typedef std::map< Key, Values > Config;
        typedef std::map< Key, size_t > Default;

        CfgParser( const char * fileName );
        ~CfgParser() {};

        /** open a config file
         */
        void open( const char * fileName );
        /** parse the config file
         */
        void parse();

        /** check if `s' is a key and in case assign it to `k'
         */
        void key(   std::string & s, Key & k );
        /** check if `s' is a value and in case assign it to `v'
         */
        void value( std::string & s, Val & v );
        /** get the values `v' for the key `k'
         */
        const Values & values( Key & k ) const;
        void values( Key & k, Values & v ) const;
        void values( const char * k, Values & v ) const;
        /** added by GHM: return true if the key `k' exists
	 */
        bool hasKey( const char* k ) const;  
        /** get the index of the default value for the key `k'
         */
        size_t defaultValue( Key & k ) const;
        size_t defaultValue( const char * k ) const;

        /** strip comments, blanks and tabs from a string s
         */
        void strip( std::string & s );
        /** find default attribute ( a `*' at the end of the value
         *  and strip it if found
         */
        bool stripDefault( std::string & s );
        /** print the parsed [key, values]
         */
        void print( std::ostream &o );

        static const std::string commentSeparator;
        static const std::string keySeparatorLeft;
        static const std::string keySeparatorRight;
        static const std::string defaultAttribute;

    private:
        std::ifstream ifile_;
        Config cfg_;
        Default default_;
        std::string currentKey_;
        bool stripped_;
        std::string fileName_;
};

#endif
