// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "./ktTrackEff.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_ktTrackEff(void *p = 0);
   static void *newArray_ktTrackEff(Long_t size, void *p);
   static void delete_ktTrackEff(void *p);
   static void deleteArray_ktTrackEff(void *p);
   static void destruct_ktTrackEff(void *p);
   static void streamer_ktTrackEff(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ktTrackEff*)
   {
      ::ktTrackEff *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ktTrackEff >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ktTrackEff", ::ktTrackEff::Class_Version(), "ktTrackEff.hh", 20,
                  typeid(::ktTrackEff), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ktTrackEff::Dictionary, isa_proxy, 16,
                  sizeof(::ktTrackEff) );
      instance.SetNew(&new_ktTrackEff);
      instance.SetNewArray(&newArray_ktTrackEff);
      instance.SetDelete(&delete_ktTrackEff);
      instance.SetDeleteArray(&deleteArray_ktTrackEff);
      instance.SetDestructor(&destruct_ktTrackEff);
      instance.SetStreamerFunc(&streamer_ktTrackEff);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ktTrackEff*)
   {
      return GenerateInitInstanceLocal((::ktTrackEff*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ktTrackEff*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ktTrackEff::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ktTrackEff::Class_Name()
{
   return "ktTrackEff";
}

//______________________________________________________________________________
const char *ktTrackEff::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ktTrackEff*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ktTrackEff::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ktTrackEff*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ktTrackEff::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ktTrackEff*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ktTrackEff::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ktTrackEff*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void ktTrackEff::Streamer(TBuffer &R__b)
{
   // Stream an object of class ktTrackEff.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      fName.Streamer(R__b);
      int R__i;
      for (R__i = 0; R__i < 3; R__i++)
         R__b >> effY04[R__i];
      for (R__i = 0; R__i < 3; R__i++)
         R__b >> effY07pteta[R__i];
      for (R__i = 0; R__i < 3; R__i++)
         R__b >> effY07eta[R__i];
      R__b >> effY06;
      R__b >> sysUn;
      R__b.CheckByteCount(R__s, R__c, ktTrackEff::IsA());
   } else {
      R__c = R__b.WriteVersion(ktTrackEff::IsA(), kTRUE);
      TObject::Streamer(R__b);
      fName.Streamer(R__b);
      int R__i;
      for (R__i = 0; R__i < 3; R__i++)
         R__b << effY04[R__i];
      for (R__i = 0; R__i < 3; R__i++)
         R__b << (TObject*)effY07pteta[R__i];
      for (R__i = 0; R__i < 3; R__i++)
         R__b << (TObject*)effY07eta[R__i];
      R__b << effY06;
      R__b << sysUn;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ktTrackEff(void *p) {
      return  p ? new(p) ::ktTrackEff : new ::ktTrackEff;
   }
   static void *newArray_ktTrackEff(Long_t nElements, void *p) {
      return p ? new(p) ::ktTrackEff[nElements] : new ::ktTrackEff[nElements];
   }
   // Wrapper around operator delete
   static void delete_ktTrackEff(void *p) {
      delete ((::ktTrackEff*)p);
   }
   static void deleteArray_ktTrackEff(void *p) {
      delete [] ((::ktTrackEff*)p);
   }
   static void destruct_ktTrackEff(void *p) {
      typedef ::ktTrackEff current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_ktTrackEff(TBuffer &buf, void *obj) {
      ((::ktTrackEff*)obj)->::ktTrackEff::Streamer(buf);
   }
} // end of namespace ROOT for class ::ktTrackEff

namespace {
  void TriggerDictionaryInitialization_dict_Impl() {
    static const char* headers[] = {
"./ktTrackEff.hh",
0
    };
    static const char* includePaths[] = {
"/opt/local/libexec/root6/include/root",
"/Users/vverkest/physics/analysis/detector_vs_particle/src/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$./ktTrackEff.hh")))  ktTrackEff;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "./ktTrackEff.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"ktTrackEff", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dict() {
  TriggerDictionaryInitialization_dict_Impl();
}
