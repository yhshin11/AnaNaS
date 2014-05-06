#include <cassert>
#include <map>
#include <string>
#include <typeinfo>
#include <cmath>
#include <sstream>
#include <iostream>
using namespace std;

#include <TFile.h>
#include <TBranch.h>
#include <TKey.h>
#include <TTree.h>
#include <TString.h>
#include <TBits.h>

#include "Analysis/utils/CutUtils.hh"

int main()
{
  cout << "hello. World" << endl;

  float var1(15);  float cut1(10);
  float var2(15);  float cut2(10); 
  float var3(15);  float cut3(10);
  float var4(15);  float cut4(10);
  float var5(15);  float cut5(10);
  float var6(15);  float cut6(10);
  
  OneSidedCut< float > c1("cut1", var1, cut1, ">" );
  OneSidedCut< float > c2("cut2", var2, cut2, "<" );
  OneSidedCut< float > c3("cut3", var3, cut3, ">" );
  OneSidedCut< float > c4("cut4", var4, cut4, ">" );
  OneSidedCut< float > c5("cut5", var5, cut5, ">" );
  OneSidedCut< float > c6("cut6", var6, cut6, ">" );
  CompositeCut cc1("cc1","&&");
  cc1.add(c1);
  cc1.add(c2);
  CompositeCut cc2("cc2","&&");
  cc2.add(cc1);
  cc2.add(c3);
  cc2.add(c4);
  CompositeCut cc3("cc3","&&");
  cc3.add(c5);
  cc3.add(c6);
  CompositeCut cc("cc","&&");
  //  cc.add(cc1);
  cc.add(cc2);
  cc.add(cc3);

  TBits bits;
  string str;
  cc.getBits(bits);
  cc.getString(str);
  cc.print();
  cout << str << ": " << bits << endl;

  cout << "N-1: cc[" << cc1.name() << "]=" << cc.ignore(cc1) << endl;
  cout << "N-1: cc[" << cc2.name() << "]=" << cc.ignore(cc2) << endl;
  cout << "N-1: cc[" << cc3.name() << "]=" << cc.ignore(cc3) << endl;
  cout << "N-1: cc[" << c1.name() << "]=" << cc.ignore(c1) << endl;
  cout << "N-1: cc[" << c2.name() << "]=" << cc.ignore(c2) << endl;
  cout << "N-1: cc[" << c3.name() << "]=" << cc.ignore(c3) << endl;
  cout << "N-1: cc[" << c4.name() << "]=" << cc.ignore(c4) << endl;
  cout << "N-1: cc[" << c5.name() << "]=" << cc.ignore(c5) << endl;
  cout << "N-1: cc[" << c6.name() << "]=" << cc.ignore(c6) << endl;

  //  bool ok = true;
  cout << "SEQ: cc[" << cc1.name() << "]=" << cc.sequential(cc1) << endl;
  cout << "SEQ: cc[" << cc2.name() << "]=" << cc.sequential(cc2) << endl;
  cout << "SEQ: cc[" << cc3.name() << "]=" << cc.sequential(cc3) << endl;
  cout << "SEQ: cc[" << c1.name() << "]=" << cc.sequential(c1) << endl;
  cout << "SEQ: cc[" << c2.name() << "]=" << cc.sequential(c2) << endl;
  cout << "SEQ: cc[" << c3.name() << "]=" << cc.sequential(c3) << endl;
  cout << "SEQ: cc[" << c4.name() << "]=" << cc.sequential(c4) << endl;
  cout << "SEQ: cc[" << c5.name() << "]=" << cc.sequential(c5) << endl;
  cout << "SEQ: cc[" << c6.name() << "]=" << cc.sequential(c6) << endl;
  

  return 0;
}
