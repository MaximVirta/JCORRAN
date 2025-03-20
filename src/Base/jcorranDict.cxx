// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME jcorranDict
#define R__NO_DEPRECATION

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

// Header files passed as explicit arguments
#include "../../include/AliAnalysisAnaTwoMultiCorrelations.h"
#include "../../include/AliAnalysisPtVn.h"
#include "../../include/AliAnalysisSPCRun2.h"
#include "../../include/AliJBaseEventHeader.h"
#include "../../include/AliJBaseTrack.h"
#include "../../include/AliJConst.h"
#include "../../include/AliJEfficiency.h"
#include "../../include/AliJFFlucAnalysis.h"
#include "../../include/AliJHistManager.h"
#include "../../include/AliJRunTable.h"
#include "../../include/AliJTrackCut.h"
#include "../../include/AliJTrack.h"
#include "../../include/JTreeDataManager.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *AliJTrackCut_Dictionary();
   static void AliJTrackCut_TClassManip(TClass*);
   static void *new_AliJTrackCut(void *p = nullptr);
   static void *newArray_AliJTrackCut(Long_t size, void *p);
   static void delete_AliJTrackCut(void *p);
   static void deleteArray_AliJTrackCut(void *p);
   static void destruct_AliJTrackCut(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJTrackCut*)
   {
      ::AliJTrackCut *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliJTrackCut));
      static ::ROOT::TGenericClassInfo 
         instance("AliJTrackCut", "../../include/AliJTrackCut.h", 11,
                  typeid(::AliJTrackCut), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliJTrackCut_Dictionary, isa_proxy, 0,
                  sizeof(::AliJTrackCut) );
      instance.SetNew(&new_AliJTrackCut);
      instance.SetNewArray(&newArray_AliJTrackCut);
      instance.SetDelete(&delete_AliJTrackCut);
      instance.SetDeleteArray(&deleteArray_AliJTrackCut);
      instance.SetDestructor(&destruct_AliJTrackCut);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJTrackCut*)
   {
      return GenerateInitInstanceLocal((::AliJTrackCut*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJTrackCut*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliJTrackCut_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliJTrackCut*)nullptr)->GetClass();
      AliJTrackCut_TClassManip(theClass);
   return theClass;
   }

   static void AliJTrackCut_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliJRunTable_Dictionary();
   static void AliJRunTable_TClassManip(TClass*);
   static void *new_AliJRunTable(void *p = nullptr);
   static void *newArray_AliJRunTable(Long_t size, void *p);
   static void delete_AliJRunTable(void *p);
   static void deleteArray_AliJRunTable(void *p);
   static void destruct_AliJRunTable(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJRunTable*)
   {
      ::AliJRunTable *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliJRunTable));
      static ::ROOT::TGenericClassInfo 
         instance("AliJRunTable", "../../include/AliJRunTable.h", 11,
                  typeid(::AliJRunTable), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliJRunTable_Dictionary, isa_proxy, 0,
                  sizeof(::AliJRunTable) );
      instance.SetNew(&new_AliJRunTable);
      instance.SetNewArray(&newArray_AliJRunTable);
      instance.SetDelete(&delete_AliJRunTable);
      instance.SetDeleteArray(&deleteArray_AliJRunTable);
      instance.SetDestructor(&destruct_AliJRunTable);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJRunTable*)
   {
      return GenerateInitInstanceLocal((::AliJRunTable*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJRunTable*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliJRunTable_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliJRunTable*)nullptr)->GetClass();
      AliJRunTable_TClassManip(theClass);
   return theClass;
   }

   static void AliJRunTable_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliJEfficiency_Dictionary();
   static void AliJEfficiency_TClassManip(TClass*);
   static void *new_AliJEfficiency(void *p = nullptr);
   static void *newArray_AliJEfficiency(Long_t size, void *p);
   static void delete_AliJEfficiency(void *p);
   static void deleteArray_AliJEfficiency(void *p);
   static void destruct_AliJEfficiency(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJEfficiency*)
   {
      ::AliJEfficiency *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliJEfficiency));
      static ::ROOT::TGenericClassInfo 
         instance("AliJEfficiency", "../../include/AliJEfficiency.h", 17,
                  typeid(::AliJEfficiency), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliJEfficiency_Dictionary, isa_proxy, 0,
                  sizeof(::AliJEfficiency) );
      instance.SetNew(&new_AliJEfficiency);
      instance.SetNewArray(&newArray_AliJEfficiency);
      instance.SetDelete(&delete_AliJEfficiency);
      instance.SetDeleteArray(&deleteArray_AliJEfficiency);
      instance.SetDestructor(&destruct_AliJEfficiency);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJEfficiency*)
   {
      return GenerateInitInstanceLocal((::AliJEfficiency*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJEfficiency*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliJEfficiency_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliJEfficiency*)nullptr)->GetClass();
      AliJEfficiency_TClassManip(theClass);
   return theClass;
   }

   static void AliJEfficiency_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliAnalysisAnaTwoMultiCorrelations_Dictionary();
   static void AliAnalysisAnaTwoMultiCorrelations_TClassManip(TClass*);
   static void *new_AliAnalysisAnaTwoMultiCorrelations(void *p = nullptr);
   static void *newArray_AliAnalysisAnaTwoMultiCorrelations(Long_t size, void *p);
   static void delete_AliAnalysisAnaTwoMultiCorrelations(void *p);
   static void deleteArray_AliAnalysisAnaTwoMultiCorrelations(void *p);
   static void destruct_AliAnalysisAnaTwoMultiCorrelations(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisAnaTwoMultiCorrelations*)
   {
      ::AliAnalysisAnaTwoMultiCorrelations *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliAnalysisAnaTwoMultiCorrelations));
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisAnaTwoMultiCorrelations", "../../include/AliAnalysisAnaTwoMultiCorrelations.h", 27,
                  typeid(::AliAnalysisAnaTwoMultiCorrelations), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliAnalysisAnaTwoMultiCorrelations_Dictionary, isa_proxy, 0,
                  sizeof(::AliAnalysisAnaTwoMultiCorrelations) );
      instance.SetNew(&new_AliAnalysisAnaTwoMultiCorrelations);
      instance.SetNewArray(&newArray_AliAnalysisAnaTwoMultiCorrelations);
      instance.SetDelete(&delete_AliAnalysisAnaTwoMultiCorrelations);
      instance.SetDeleteArray(&deleteArray_AliAnalysisAnaTwoMultiCorrelations);
      instance.SetDestructor(&destruct_AliAnalysisAnaTwoMultiCorrelations);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisAnaTwoMultiCorrelations*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisAnaTwoMultiCorrelations*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisAnaTwoMultiCorrelations*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliAnalysisAnaTwoMultiCorrelations_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisAnaTwoMultiCorrelations*)nullptr)->GetClass();
      AliAnalysisAnaTwoMultiCorrelations_TClassManip(theClass);
   return theClass;
   }

   static void AliAnalysisAnaTwoMultiCorrelations_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_AliJBaseTrack(void *p = nullptr);
   static void *newArray_AliJBaseTrack(Long_t size, void *p);
   static void delete_AliJBaseTrack(void *p);
   static void deleteArray_AliJBaseTrack(void *p);
   static void destruct_AliJBaseTrack(void *p);
   static void streamer_AliJBaseTrack(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJBaseTrack*)
   {
      ::AliJBaseTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::AliJBaseTrack >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("AliJBaseTrack", ::AliJBaseTrack::Class_Version(), "../../include/AliJBaseTrack.h", 33,
                  typeid(::AliJBaseTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::AliJBaseTrack::Dictionary, isa_proxy, 16,
                  sizeof(::AliJBaseTrack) );
      instance.SetNew(&new_AliJBaseTrack);
      instance.SetNewArray(&newArray_AliJBaseTrack);
      instance.SetDelete(&delete_AliJBaseTrack);
      instance.SetDeleteArray(&deleteArray_AliJBaseTrack);
      instance.SetDestructor(&destruct_AliJBaseTrack);
      instance.SetStreamerFunc(&streamer_AliJBaseTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJBaseTrack*)
   {
      return GenerateInitInstanceLocal((::AliJBaseTrack*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJBaseTrack*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *AliAnalysisPtVn_Dictionary();
   static void AliAnalysisPtVn_TClassManip(TClass*);
   static void *new_AliAnalysisPtVn(void *p = nullptr);
   static void *newArray_AliAnalysisPtVn(Long_t size, void *p);
   static void delete_AliAnalysisPtVn(void *p);
   static void deleteArray_AliAnalysisPtVn(void *p);
   static void destruct_AliAnalysisPtVn(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisPtVn*)
   {
      ::AliAnalysisPtVn *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliAnalysisPtVn));
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisPtVn", "../../include/AliAnalysisPtVn.h", 29,
                  typeid(::AliAnalysisPtVn), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliAnalysisPtVn_Dictionary, isa_proxy, 0,
                  sizeof(::AliAnalysisPtVn) );
      instance.SetNew(&new_AliAnalysisPtVn);
      instance.SetNewArray(&newArray_AliAnalysisPtVn);
      instance.SetDelete(&delete_AliAnalysisPtVn);
      instance.SetDeleteArray(&deleteArray_AliAnalysisPtVn);
      instance.SetDestructor(&destruct_AliAnalysisPtVn);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisPtVn*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisPtVn*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisPtVn*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliAnalysisPtVn_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisPtVn*)nullptr)->GetClass();
      AliAnalysisPtVn_TClassManip(theClass);
   return theClass;
   }

   static void AliAnalysisPtVn_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliAnalysisSPCRun2_Dictionary();
   static void AliAnalysisSPCRun2_TClassManip(TClass*);
   static void *new_AliAnalysisSPCRun2(void *p = nullptr);
   static void *newArray_AliAnalysisSPCRun2(Long_t size, void *p);
   static void delete_AliAnalysisSPCRun2(void *p);
   static void deleteArray_AliAnalysisSPCRun2(void *p);
   static void destruct_AliAnalysisSPCRun2(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisSPCRun2*)
   {
      ::AliAnalysisSPCRun2 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliAnalysisSPCRun2));
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisSPCRun2", "../../include/AliAnalysisSPCRun2.h", 27,
                  typeid(::AliAnalysisSPCRun2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliAnalysisSPCRun2_Dictionary, isa_proxy, 0,
                  sizeof(::AliAnalysisSPCRun2) );
      instance.SetNew(&new_AliAnalysisSPCRun2);
      instance.SetNewArray(&newArray_AliAnalysisSPCRun2);
      instance.SetDelete(&delete_AliAnalysisSPCRun2);
      instance.SetDeleteArray(&deleteArray_AliAnalysisSPCRun2);
      instance.SetDestructor(&destruct_AliAnalysisSPCRun2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisSPCRun2*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisSPCRun2*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisSPCRun2*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliAnalysisSPCRun2_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisSPCRun2*)nullptr)->GetClass();
      AliAnalysisSPCRun2_TClassManip(theClass);
   return theClass;
   }

   static void AliAnalysisSPCRun2_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_AliJBaseEventHeader(void *p = nullptr);
   static void *newArray_AliJBaseEventHeader(Long_t size, void *p);
   static void delete_AliJBaseEventHeader(void *p);
   static void deleteArray_AliJBaseEventHeader(void *p);
   static void destruct_AliJBaseEventHeader(void *p);
   static void streamer_AliJBaseEventHeader(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJBaseEventHeader*)
   {
      ::AliJBaseEventHeader *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::AliJBaseEventHeader >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("AliJBaseEventHeader", ::AliJBaseEventHeader::Class_Version(), "../../include/AliJBaseEventHeader.h", 27,
                  typeid(::AliJBaseEventHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::AliJBaseEventHeader::Dictionary, isa_proxy, 16,
                  sizeof(::AliJBaseEventHeader) );
      instance.SetNew(&new_AliJBaseEventHeader);
      instance.SetNewArray(&newArray_AliJBaseEventHeader);
      instance.SetDelete(&delete_AliJBaseEventHeader);
      instance.SetDeleteArray(&deleteArray_AliJBaseEventHeader);
      instance.SetDestructor(&destruct_AliJBaseEventHeader);
      instance.SetStreamerFunc(&streamer_AliJBaseEventHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJBaseEventHeader*)
   {
      return GenerateInitInstanceLocal((::AliJBaseEventHeader*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJBaseEventHeader*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *AliJHistManager_Dictionary();
   static void AliJHistManager_TClassManip(TClass*);
   static void delete_AliJHistManager(void *p);
   static void deleteArray_AliJHistManager(void *p);
   static void destruct_AliJHistManager(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJHistManager*)
   {
      ::AliJHistManager *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliJHistManager));
      static ::ROOT::TGenericClassInfo 
         instance("AliJHistManager", "../../include/AliJHistManager.h", 386,
                  typeid(::AliJHistManager), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliJHistManager_Dictionary, isa_proxy, 0,
                  sizeof(::AliJHistManager) );
      instance.SetDelete(&delete_AliJHistManager);
      instance.SetDeleteArray(&deleteArray_AliJHistManager);
      instance.SetDestructor(&destruct_AliJHistManager);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJHistManager*)
   {
      return GenerateInitInstanceLocal((::AliJHistManager*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJHistManager*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliJHistManager_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliJHistManager*)nullptr)->GetClass();
      AliJHistManager_TClassManip(theClass);
   return theClass;
   }

   static void AliJHistManager_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliJFFlucAnalysis_Dictionary();
   static void AliJFFlucAnalysis_TClassManip(TClass*);
   static void *new_AliJFFlucAnalysis(void *p = nullptr);
   static void *newArray_AliJFFlucAnalysis(Long_t size, void *p);
   static void delete_AliJFFlucAnalysis(void *p);
   static void deleteArray_AliJFFlucAnalysis(void *p);
   static void destruct_AliJFFlucAnalysis(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJFFlucAnalysis*)
   {
      ::AliJFFlucAnalysis *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliJFFlucAnalysis));
      static ::ROOT::TGenericClassInfo 
         instance("AliJFFlucAnalysis", "../../include/AliJFFlucAnalysis.h", 13,
                  typeid(::AliJFFlucAnalysis), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliJFFlucAnalysis_Dictionary, isa_proxy, 0,
                  sizeof(::AliJFFlucAnalysis) );
      instance.SetNew(&new_AliJFFlucAnalysis);
      instance.SetNewArray(&newArray_AliJFFlucAnalysis);
      instance.SetDelete(&delete_AliJFFlucAnalysis);
      instance.SetDeleteArray(&deleteArray_AliJFFlucAnalysis);
      instance.SetDestructor(&destruct_AliJFFlucAnalysis);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJFFlucAnalysis*)
   {
      return GenerateInitInstanceLocal((::AliJFFlucAnalysis*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJFFlucAnalysis*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliJFFlucAnalysis_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliJFFlucAnalysis*)nullptr)->GetClass();
      AliJFFlucAnalysis_TClassManip(theClass);
   return theClass;
   }

   static void AliJFFlucAnalysis_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_AliJTrack(void *p = nullptr);
   static void *newArray_AliJTrack(Long_t size, void *p);
   static void delete_AliJTrack(void *p);
   static void deleteArray_AliJTrack(void *p);
   static void destruct_AliJTrack(void *p);
   static void streamer_AliJTrack(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliJTrack*)
   {
      ::AliJTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::AliJTrack >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("AliJTrack", ::AliJTrack::Class_Version(), "../../include/AliJTrack.h", 27,
                  typeid(::AliJTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::AliJTrack::Dictionary, isa_proxy, 16,
                  sizeof(::AliJTrack) );
      instance.SetNew(&new_AliJTrack);
      instance.SetNewArray(&newArray_AliJTrack);
      instance.SetDelete(&delete_AliJTrack);
      instance.SetDeleteArray(&deleteArray_AliJTrack);
      instance.SetDestructor(&destruct_AliJTrack);
      instance.SetStreamerFunc(&streamer_AliJTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliJTrack*)
   {
      return GenerateInitInstanceLocal((::AliJTrack*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliJTrack*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *JTreeDataManager_Dictionary();
   static void JTreeDataManager_TClassManip(TClass*);
   static void *new_JTreeDataManager(void *p = nullptr);
   static void *newArray_JTreeDataManager(Long_t size, void *p);
   static void delete_JTreeDataManager(void *p);
   static void deleteArray_JTreeDataManager(void *p);
   static void destruct_JTreeDataManager(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JTreeDataManager*)
   {
      ::JTreeDataManager *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::JTreeDataManager));
      static ::ROOT::TGenericClassInfo 
         instance("JTreeDataManager", "../../include/JTreeDataManager.h", 36,
                  typeid(::JTreeDataManager), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &JTreeDataManager_Dictionary, isa_proxy, 0,
                  sizeof(::JTreeDataManager) );
      instance.SetNew(&new_JTreeDataManager);
      instance.SetNewArray(&newArray_JTreeDataManager);
      instance.SetDelete(&delete_JTreeDataManager);
      instance.SetDeleteArray(&deleteArray_JTreeDataManager);
      instance.SetDestructor(&destruct_JTreeDataManager);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JTreeDataManager*)
   {
      return GenerateInitInstanceLocal((::JTreeDataManager*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::JTreeDataManager*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *JTreeDataManager_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::JTreeDataManager*)nullptr)->GetClass();
      JTreeDataManager_TClassManip(theClass);
   return theClass;
   }

   static void JTreeDataManager_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr AliJBaseTrack::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *AliJBaseTrack::Class_Name()
{
   return "AliJBaseTrack";
}

//______________________________________________________________________________
const char *AliJBaseTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliJBaseTrack*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int AliJBaseTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliJBaseTrack*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *AliJBaseTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliJBaseTrack*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *AliJBaseTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliJBaseTrack*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr AliJBaseEventHeader::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *AliJBaseEventHeader::Class_Name()
{
   return "AliJBaseEventHeader";
}

//______________________________________________________________________________
const char *AliJBaseEventHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliJBaseEventHeader*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int AliJBaseEventHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliJBaseEventHeader*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *AliJBaseEventHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliJBaseEventHeader*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *AliJBaseEventHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliJBaseEventHeader*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr AliJTrack::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *AliJTrack::Class_Name()
{
   return "AliJTrack";
}

//______________________________________________________________________________
const char *AliJTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliJTrack*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int AliJTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliJTrack*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *AliJTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliJTrack*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *AliJTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliJTrack*)nullptr)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJTrackCut(void *p) {
      return  p ? new(p) ::AliJTrackCut : new ::AliJTrackCut;
   }
   static void *newArray_AliJTrackCut(Long_t nElements, void *p) {
      return p ? new(p) ::AliJTrackCut[nElements] : new ::AliJTrackCut[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJTrackCut(void *p) {
      delete ((::AliJTrackCut*)p);
   }
   static void deleteArray_AliJTrackCut(void *p) {
      delete [] ((::AliJTrackCut*)p);
   }
   static void destruct_AliJTrackCut(void *p) {
      typedef ::AliJTrackCut current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliJTrackCut

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJRunTable(void *p) {
      return  p ? new(p) ::AliJRunTable : new ::AliJRunTable;
   }
   static void *newArray_AliJRunTable(Long_t nElements, void *p) {
      return p ? new(p) ::AliJRunTable[nElements] : new ::AliJRunTable[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJRunTable(void *p) {
      delete ((::AliJRunTable*)p);
   }
   static void deleteArray_AliJRunTable(void *p) {
      delete [] ((::AliJRunTable*)p);
   }
   static void destruct_AliJRunTable(void *p) {
      typedef ::AliJRunTable current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliJRunTable

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJEfficiency(void *p) {
      return  p ? new(p) ::AliJEfficiency : new ::AliJEfficiency;
   }
   static void *newArray_AliJEfficiency(Long_t nElements, void *p) {
      return p ? new(p) ::AliJEfficiency[nElements] : new ::AliJEfficiency[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJEfficiency(void *p) {
      delete ((::AliJEfficiency*)p);
   }
   static void deleteArray_AliJEfficiency(void *p) {
      delete [] ((::AliJEfficiency*)p);
   }
   static void destruct_AliJEfficiency(void *p) {
      typedef ::AliJEfficiency current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliJEfficiency

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisAnaTwoMultiCorrelations(void *p) {
      return  p ? new(p) ::AliAnalysisAnaTwoMultiCorrelations : new ::AliAnalysisAnaTwoMultiCorrelations;
   }
   static void *newArray_AliAnalysisAnaTwoMultiCorrelations(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisAnaTwoMultiCorrelations[nElements] : new ::AliAnalysisAnaTwoMultiCorrelations[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisAnaTwoMultiCorrelations(void *p) {
      delete ((::AliAnalysisAnaTwoMultiCorrelations*)p);
   }
   static void deleteArray_AliAnalysisAnaTwoMultiCorrelations(void *p) {
      delete [] ((::AliAnalysisAnaTwoMultiCorrelations*)p);
   }
   static void destruct_AliAnalysisAnaTwoMultiCorrelations(void *p) {
      typedef ::AliAnalysisAnaTwoMultiCorrelations current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisAnaTwoMultiCorrelations

//______________________________________________________________________________
void AliJBaseTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliJBaseTrack.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TLorentzVector::Streamer(R__b);
      R__b >> fID;
      R__b >> fLabel;
      R__b >> fParticleType;
      R__b >> fCharge;
      R__b >> fStatus;
      R__b >> fFlags;
      R__b.CheckByteCount(R__s, R__c, AliJBaseTrack::IsA());
   } else {
      R__c = R__b.WriteVersion(AliJBaseTrack::IsA(), kTRUE);
      TLorentzVector::Streamer(R__b);
      R__b << fID;
      R__b << fLabel;
      R__b << fParticleType;
      R__b << fCharge;
      R__b << fStatus;
      R__b << fFlags;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJBaseTrack(void *p) {
      return  p ? new(p) ::AliJBaseTrack : new ::AliJBaseTrack;
   }
   static void *newArray_AliJBaseTrack(Long_t nElements, void *p) {
      return p ? new(p) ::AliJBaseTrack[nElements] : new ::AliJBaseTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJBaseTrack(void *p) {
      delete ((::AliJBaseTrack*)p);
   }
   static void deleteArray_AliJBaseTrack(void *p) {
      delete [] ((::AliJBaseTrack*)p);
   }
   static void destruct_AliJBaseTrack(void *p) {
      typedef ::AliJBaseTrack current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_AliJBaseTrack(TBuffer &buf, void *obj) {
      ((::AliJBaseTrack*)obj)->::AliJBaseTrack::Streamer(buf);
   }
} // end of namespace ROOT for class ::AliJBaseTrack

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisPtVn(void *p) {
      return  p ? new(p) ::AliAnalysisPtVn : new ::AliAnalysisPtVn;
   }
   static void *newArray_AliAnalysisPtVn(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisPtVn[nElements] : new ::AliAnalysisPtVn[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisPtVn(void *p) {
      delete ((::AliAnalysisPtVn*)p);
   }
   static void deleteArray_AliAnalysisPtVn(void *p) {
      delete [] ((::AliAnalysisPtVn*)p);
   }
   static void destruct_AliAnalysisPtVn(void *p) {
      typedef ::AliAnalysisPtVn current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisPtVn

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisSPCRun2(void *p) {
      return  p ? new(p) ::AliAnalysisSPCRun2 : new ::AliAnalysisSPCRun2;
   }
   static void *newArray_AliAnalysisSPCRun2(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisSPCRun2[nElements] : new ::AliAnalysisSPCRun2[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisSPCRun2(void *p) {
      delete ((::AliAnalysisSPCRun2*)p);
   }
   static void deleteArray_AliAnalysisSPCRun2(void *p) {
      delete [] ((::AliAnalysisSPCRun2*)p);
   }
   static void destruct_AliAnalysisSPCRun2(void *p) {
      typedef ::AliAnalysisSPCRun2 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisSPCRun2

//______________________________________________________________________________
void AliJBaseEventHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliJBaseEventHeader.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> fEventID;
      {float R_Dummy; R__b >> R_Dummy; fCentrality=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxX=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxY=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxZ=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxZErr=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxMCX=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxMCY=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxMCZ=Double32_t(R_Dummy);}
      R__b.CheckByteCount(R__s, R__c, AliJBaseEventHeader::IsA());
   } else {
      R__c = R__b.WriteVersion(AliJBaseEventHeader::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << fEventID;
      R__b << float(fCentrality);
      R__b << float(fVtxX);
      R__b << float(fVtxY);
      R__b << float(fVtxZ);
      R__b << float(fVtxZErr);
      R__b << float(fVtxMCX);
      R__b << float(fVtxMCY);
      R__b << float(fVtxMCZ);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJBaseEventHeader(void *p) {
      return  p ? new(p) ::AliJBaseEventHeader : new ::AliJBaseEventHeader;
   }
   static void *newArray_AliJBaseEventHeader(Long_t nElements, void *p) {
      return p ? new(p) ::AliJBaseEventHeader[nElements] : new ::AliJBaseEventHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJBaseEventHeader(void *p) {
      delete ((::AliJBaseEventHeader*)p);
   }
   static void deleteArray_AliJBaseEventHeader(void *p) {
      delete [] ((::AliJBaseEventHeader*)p);
   }
   static void destruct_AliJBaseEventHeader(void *p) {
      typedef ::AliJBaseEventHeader current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_AliJBaseEventHeader(TBuffer &buf, void *obj) {
      ((::AliJBaseEventHeader*)obj)->::AliJBaseEventHeader::Streamer(buf);
   }
} // end of namespace ROOT for class ::AliJBaseEventHeader

namespace ROOT {
   // Wrapper around operator delete
   static void delete_AliJHistManager(void *p) {
      delete ((::AliJHistManager*)p);
   }
   static void deleteArray_AliJHistManager(void *p) {
      delete [] ((::AliJHistManager*)p);
   }
   static void destruct_AliJHistManager(void *p) {
      typedef ::AliJHistManager current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliJHistManager

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJFFlucAnalysis(void *p) {
      return  p ? new(p) ::AliJFFlucAnalysis : new ::AliJFFlucAnalysis;
   }
   static void *newArray_AliJFFlucAnalysis(Long_t nElements, void *p) {
      return p ? new(p) ::AliJFFlucAnalysis[nElements] : new ::AliJFFlucAnalysis[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJFFlucAnalysis(void *p) {
      delete ((::AliJFFlucAnalysis*)p);
   }
   static void deleteArray_AliJFFlucAnalysis(void *p) {
      delete [] ((::AliJFFlucAnalysis*)p);
   }
   static void destruct_AliJFFlucAnalysis(void *p) {
      typedef ::AliJFFlucAnalysis current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliJFFlucAnalysis

//______________________________________________________________________________
void AliJTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliJTrack.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      AliJBaseTrack::Streamer(R__b);
      R__b.ReadStaticArrayDouble32((double*)fTrkPID);
      R__b >> fFilterMap;
      R__b >> fTPCnClust;
      {float R_Dummy; R__b >> R_Dummy; fTPCdEdx=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fTOFsignal=Double32_t(R_Dummy);}
      R__b.ReadStaticArrayDouble32(fExpTOFsignal);
      {float R_Dummy; R__b >> R_Dummy; fTPCmom=Double32_t(R_Dummy);}
      R__b.ReadStaticArrayDouble32(fTPCTrack);
      R__b.ReadStaticArrayDouble32(fGCGTrack);
      R__b.ReadStaticArrayDouble32(fTrackPos);
      R__b.CheckByteCount(R__s, R__c, AliJTrack::IsA());
   } else {
      R__c = R__b.WriteVersion(AliJTrack::IsA(), kTRUE);
      AliJBaseTrack::Streamer(R__b);
      R__b.WriteArrayDouble32((double*)fTrkPID, 5);
      R__b << fFilterMap;
      R__b << fTPCnClust;
      R__b << float(fTPCdEdx);
      R__b << float(fTOFsignal);
      R__b.WriteArrayDouble32(fExpTOFsignal, 5);
      R__b << float(fTPCmom);
      R__b.WriteArrayDouble32(fTPCTrack, 3);
      R__b.WriteArrayDouble32(fGCGTrack, 3);
      R__b.WriteArrayDouble32(fTrackPos, 3);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliJTrack(void *p) {
      return  p ? new(p) ::AliJTrack : new ::AliJTrack;
   }
   static void *newArray_AliJTrack(Long_t nElements, void *p) {
      return p ? new(p) ::AliJTrack[nElements] : new ::AliJTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliJTrack(void *p) {
      delete ((::AliJTrack*)p);
   }
   static void deleteArray_AliJTrack(void *p) {
      delete [] ((::AliJTrack*)p);
   }
   static void destruct_AliJTrack(void *p) {
      typedef ::AliJTrack current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_AliJTrack(TBuffer &buf, void *obj) {
      ((::AliJTrack*)obj)->::AliJTrack::Streamer(buf);
   }
} // end of namespace ROOT for class ::AliJTrack

namespace ROOT {
   // Wrappers around operator new
   static void *new_JTreeDataManager(void *p) {
      return  p ? new(p) ::JTreeDataManager : new ::JTreeDataManager;
   }
   static void *newArray_JTreeDataManager(Long_t nElements, void *p) {
      return p ? new(p) ::JTreeDataManager[nElements] : new ::JTreeDataManager[nElements];
   }
   // Wrapper around operator delete
   static void delete_JTreeDataManager(void *p) {
      delete ((::JTreeDataManager*)p);
   }
   static void deleteArray_JTreeDataManager(void *p) {
      delete [] ((::JTreeDataManager*)p);
   }
   static void destruct_JTreeDataManager(void *p) {
      typedef ::JTreeDataManager current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::JTreeDataManager

namespace {
  void TriggerDictionaryInitialization_jcorranDict_Impl() {
    static const char* headers[] = {
"../../include/AliAnalysisAnaTwoMultiCorrelations.h",
"../../include/AliAnalysisPtVn.h",
"../../include/AliAnalysisSPCRun2.h",
"../../include/AliJBaseEventHeader.h",
"../../include/AliJBaseTrack.h",
"../../include/AliJConst.h",
"../../include/AliJEfficiency.h",
"../../include/AliJFFlucAnalysis.h",
"../../include/AliJHistManager.h",
"../../include/AliJRunTable.h",
"../../include/AliJTrackCut.h",
"../../include/AliJTrack.h",
"../../include/JTreeDataManager.h",
nullptr
    };
    static const char* includePaths[] = {
"/projappl/project_2003112/root/include/",
"/projappl/project_2003154/lib/JCORRAN/src/Base/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "jcorranDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$../../include/AliAnalysisAnaTwoMultiCorrelations.h")))  AliJTrackCut;
class __attribute__((annotate("$clingAutoload$../../include/AliAnalysisAnaTwoMultiCorrelations.h")))  AliJRunTable;
class __attribute__((annotate("$clingAutoload$../../include/AliAnalysisAnaTwoMultiCorrelations.h")))  AliJEfficiency;
class __attribute__((annotate("$clingAutoload$../../include/AliAnalysisAnaTwoMultiCorrelations.h")))  AliAnalysisAnaTwoMultiCorrelations;
class __attribute__((annotate("$clingAutoload$../../include/AliAnalysisPtVn.h")))  AliJBaseTrack;
class __attribute__((annotate("$clingAutoload$../../include/AliAnalysisPtVn.h")))  AliAnalysisPtVn;
class __attribute__((annotate("$clingAutoload$../../include/AliAnalysisSPCRun2.h")))  AliAnalysisSPCRun2;
class __attribute__((annotate("$clingAutoload$../../include/AliJBaseEventHeader.h")))  AliJBaseEventHeader;
class __attribute__((annotate("$clingAutoload$../../include/AliJFFlucAnalysis.h")))  AliJHistManager;
class __attribute__((annotate("$clingAutoload$../../include/AliJFFlucAnalysis.h")))  AliJFFlucAnalysis;
class __attribute__((annotate("$clingAutoload$../../include/AliJTrack.h")))  AliJTrack;
class __attribute__((annotate("$clingAutoload$../../include/JTreeDataManager.h")))  JTreeDataManager;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "jcorranDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "../../include/AliAnalysisAnaTwoMultiCorrelations.h"
#include "../../include/AliAnalysisPtVn.h"
#include "../../include/AliAnalysisSPCRun2.h"
#include "../../include/AliJBaseEventHeader.h"
#include "../../include/AliJBaseTrack.h"
#include "../../include/AliJConst.h"
#include "../../include/AliJEfficiency.h"
#include "../../include/AliJFFlucAnalysis.h"
#include "../../include/AliJHistManager.h"
#include "../../include/AliJRunTable.h"
#include "../../include/AliJTrackCut.h"
#include "../../include/AliJTrack.h"
#include "../../include/JTreeDataManager.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"AliAnalysisAnaTwoMultiCorrelations", payloadCode, "@",
"AliAnalysisPtVn", payloadCode, "@",
"AliAnalysisSPCRun2", payloadCode, "@",
"AliJBaseEventHeader", payloadCode, "@",
"AliJBaseTrack", payloadCode, "@",
"AliJEfficiency", payloadCode, "@",
"AliJFFlucAnalysis", payloadCode, "@",
"AliJHistManager", payloadCode, "@",
"AliJRunTable", payloadCode, "@",
"AliJTrack", payloadCode, "@",
"AliJTrackCut", payloadCode, "@",
"JTreeDataManager", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("jcorranDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_jcorranDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_jcorranDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_jcorranDict() {
  TriggerDictionaryInitialization_jcorranDict_Impl();
}
