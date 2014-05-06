/********************************************************************
* src/MECore_dict_src.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error src/MECore_dict_src.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtableMECore_dict_src();
extern void G__cpp_setup_inheritanceMECore_dict_src();
extern void G__cpp_setup_typetableMECore_dict_src();
extern void G__cpp_setup_memvarMECore_dict_src();
extern void G__cpp_setup_globalMECore_dict_src();
extern void G__cpp_setup_memfuncMECore_dict_src();
extern void G__cpp_setup_funcMECore_dict_src();
extern void G__set_cpp_environmentMECore_dict_src();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "src/ME.hh"
#include "src/MEGeom.hh"
#include "src/MEEBGeom.hh"
#include "src/MEEEGeom.hh"
#include "src/MEChannel.hh"
#include "src/MEEBDisplay.hh"
#include "src/MEEEDisplay.hh"
#include "src/MECanvasHolder.hh"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__MECore_dict_srcLN_TClass;
extern G__linked_taginfo G__MECore_dict_srcLN_TBuffer;
extern G__linked_taginfo G__MECore_dict_srcLN_TMemberInspector;
extern G__linked_taginfo G__MECore_dict_srcLN_TObject;
extern G__linked_taginfo G__MECore_dict_srcLN_TString;
extern G__linked_taginfo G__MECore_dict_srcLN_vectorlEunsignedsPlongcOallocatorlEunsignedsPlonggRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_basic_ostreamlEcharcOchar_traitslEchargRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__MECore_dict_srcLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__MECore_dict_srcLN_pairlEintcOintgR;
extern G__linked_taginfo G__MECore_dict_srcLN_MEChannel;
extern G__linked_taginfo G__MECore_dict_srcLN_ME;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLEcalRegion;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLEcalUnit;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLEcalElec;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLPN;
extern G__linked_taginfo G__MECore_dict_srcLN_vectorlEintcOallocatorlEintgRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLRunType;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLColor;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLGain;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLPNGain;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLHeader;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLSettings;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLdA;
extern G__linked_taginfo G__MECore_dict_srcLN_MEcLcLTUnit;
extern G__linked_taginfo G__MECore_dict_srcLN_vectorlEMEChannelmUcOallocatorlEMEChannelmUgRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_reverse_iteratorlEvectorlEMEChannelmUcOallocatorlEMEChannelmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__MECore_dict_srcLN_TPolyLine;
extern G__linked_taginfo G__MECore_dict_srcLN_TVectorTlEfloatgR;
extern G__linked_taginfo G__MECore_dict_srcLN_TVectorTlEdoublegR;
extern G__linked_taginfo G__MECore_dict_srcLN_TH1;
extern G__linked_taginfo G__MECore_dict_srcLN_TGraph;
extern G__linked_taginfo G__MECore_dict_srcLN_MEEBGeom;
extern G__linked_taginfo G__MECore_dict_srcLN_MEEBGeomcLcLEBUnit;
extern G__linked_taginfo G__MECore_dict_srcLN_pairlEfloatcOfloatgR;
extern G__linked_taginfo G__MECore_dict_srcLN_MEEEGeom;
extern G__linked_taginfo G__MECore_dict_srcLN_MEEEGeomcLcLEEUnit;
extern G__linked_taginfo G__MECore_dict_srcLN_listlEpairlEfloatcOfloatgRcOallocatorlEpairlEfloatcOfloatgRsPgRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_TMatrixTBaselEfloatgR;
extern G__linked_taginfo G__MECore_dict_srcLN_TMatrixTBaselEdoublegR;
extern G__linked_taginfo G__MECore_dict_srcLN_TH2;
extern G__linked_taginfo G__MECore_dict_srcLN_TCanvas;
extern G__linked_taginfo G__MECore_dict_srcLN_TPad;
extern G__linked_taginfo G__MECore_dict_srcLN_MEGeom;
extern G__linked_taginfo G__MECore_dict_srcLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR;
extern G__linked_taginfo G__MECore_dict_srcLN_MEEBDisplay;
extern G__linked_taginfo G__MECore_dict_srcLN_maplEintcOpairlEfloatcOfloatgRcOlesslEintgRcOallocatorlEpairlEconstsPintcOpairlEfloatcOfloatgRsPgRsPgRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_maplEintcOTPolyLinemUcOlesslEintgRcOallocatorlEpairlEconstsPintcOTPolyLinemUgRsPgRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_listlETObjectmUcOallocatorlETObjectmUgRsPgR;
extern G__linked_taginfo G__MECore_dict_srcLN_MEEEDisplay;
extern G__linked_taginfo G__MECore_dict_srcLN_TText;
extern G__linked_taginfo G__MECore_dict_srcLN_TLatex;
extern G__linked_taginfo G__MECore_dict_srcLN_TPaveText;
extern G__linked_taginfo G__MECore_dict_srcLN_MECanvasHolder;

/* STUB derived class for protected member access */
