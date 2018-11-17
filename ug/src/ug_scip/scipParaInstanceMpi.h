/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*             This file is part of the program and software framework       */
/*                  UG --- Ubquity Generator Framework                       */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  UG is distributed under the terms of the ZIB Academic Licence.           */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with UG; see the file COPYING. If not email to scip@zib.de.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    scipParaInstanceMpi.h
 * @brief   ScipParaInstance extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_INSTANCE_MPI_H__
#define __SCIP_PARA_INSTANCE_MPI_H__

#include <cstring>
#include <mpi.h>
#include "ug/paraDef.h"
#include "ug/paraComm.h"
#include "scipParaInstance.h"

namespace ParaSCIP
{

/** ScipInstanceMpi */
class ScipParaInstanceMpi : public ScipParaInstance
{

   int dummyToKeepStartPos;
   const char *fileName;

   /** create scipDiffParamSetPreType */
   MPI_Datatype createDatatype1();
   /** create scipDiffParamSetPreType */
   MPI_Datatype createDatatype2(bool memAllocNecessary);
   /** create scipDiffParamSetType */
   MPI_Datatype createDatatype3(bool memAllocNecessary);

   void allocateMemoryForDatatype2();
   void allocateMemoryForDatatype3();

   const char *getFileName(){ return fileName; }

public:
   /** constructor */
   ScipParaInstanceMpi(
         )
         : dummyToKeepStartPos(0), fileName(0)
   {
   }

   /** constructor : only called from ScipInitiator */
   ScipParaInstanceMpi(
         SCIP *scip,
         int  method
         ) : ScipParaInstance(scip, method), dummyToKeepStartPos(0), fileName(0)
   {
      if( method == 2 && SCIPgetStage(scip) != SCIP_STAGE_SOLVED )
      {
         SCIP *tempScip;
         SCIP_Bool success = TRUE;
         SCIP_CALL_ABORT( SCIPcreate(&tempScip) );

         /* copy all plugins and settings */
   #if SCIP_VERSION == 211 && SCIP_SUBVERSION == 0
         SCIP_CALL_ABORT( SCIPcopyPlugins(scip, tempScip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
               TRUE, TRUE, TRUE, TRUE, &success) );
   #else
      #if SCIP_APIVERSION >= 17
         SCIP_CALL_ABORT( SCIPcopyPlugins(scip, tempScip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, &success) );
      #else
         SCIP_CALL_ABORT( SCIPcopyPlugins(scip, tempScip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
               TRUE, TRUE, TRUE, TRUE, FALSE, &success) );
      #endif
   #endif
         SCIP_CALL_ABORT( SCIPcopyParamSettings(scip, tempScip) );

         /* create the variable mapping hash map */
         SCIP_HASHMAP* varmap;
         if( SCIPgetNVars(scip) > 0 )
         {  
            SCIP_CALL_ABORT( SCIPhashmapCreate(&varmap, SCIPblkmem(tempScip), SCIPgetNVars(scip)) );
         }
         SCIP_HASHMAP* conssmap;
         if( SCIPgetNConss(scip) > 0 )
         {
            SCIP_CALL_ABORT( SCIPhashmapCreate(&conssmap, SCIPblkmem(tempScip), SCIPgetNConss(scip)) );
         } 

         /* create problem in the target SCIP */
         SCIP_CALL_ABORT( SCIPcopyProb(scip, tempScip, varmap, conssmap, TRUE, "") );

         // commPth->lockApp();
         /* copy all variables and constraints */
         if( SCIPgetNVars(scip) > 0 )
         {
#if (SCIP_VERSION < 321 || ( SCIP_VERSION == 321 && SCIP_SUBVERSION < 2) )
            SCIP_CALL_ABORT( SCIPcopyVars(scip, tempScip, varmap, conssmap, TRUE) );
#else
            SCIP_CALL_ABORT( SCIPcopyVars(scip, tempScip, varmap, conssmap, NULL, NULL, 0, TRUE) );
#endif
         }
         if( SCIPgetNConss(scip) > 0 )
         {
            SCIP_CALL_ABORT( SCIPcopyConss(scip, tempScip, varmap, conssmap, TRUE, FALSE, &success) );
         }
         if( !success )
         {
            if( SCIPgetNConss(scip) > 0 )
            {
               SCIPhashmapFree(&conssmap);
            }
            if( SCIPgetNVars(scip) > 0 )
            {
               SCIPhashmapFree(&varmap);
            }
            SCIPfree(&tempScip);
            std::cerr << "Some constraint handler did not perform a valid copy. Cannot solve this instance." << std::endl;
            exit(1);
         }

         nVars = SCIPgetNVars(scip);                // original number
         int n = SCIPgetNVars(tempScip);   // copy may increase the number

         if( nVars == n )
         {
             varIndexRange = nVars;
            if( SCIPgetNConss(scip) > 0 )
            {
               SCIPhashmapFree(&conssmap);
            }
            if( SCIPgetNVars(scip) > 0 )
            {
               SCIPhashmapFree(&varmap);
            }
            SCIPfree(&tempScip);
            paraInstanceScip = scip;
            nCopies = 1;
            std::cout << "** ParaScipInstance is copied once. **" << std::endl;
         }
         else
         {
            assert(n > nVars);
            mapToOriginalIndecies = new int[SCIPgetNTotalVars(tempScip)];   // need to allocate enough for SCIPvarGetIndex(copyvar)
            mapToSolverLocalIndecies = new int[SCIPgetNTotalVars(tempScip)];
            for( int i = 0; i < SCIPgetNTotalVars(tempScip); i++ )
            {
               mapToOriginalIndecies[i] = -1;
               mapToSolverLocalIndecies[i] = -1;
            }
            SCIP_VAR **tempVars = SCIPgetVars(tempScip);
            for( int i = 0; i < n; i++ )
            {
               mapToOriginalIndecies[SCIPvarGetIndex(tempVars[i])] = i;
               mapToSolverLocalIndecies[i] = SCIPvarGetIndex(tempVars[i]);
            }

            orgScip = scip;
            nVars = n;
            varIndexRange = SCIPgetNTotalVars(tempScip);
            paraInstanceScip = tempScip;

            if( SCIPgetNConss(scip) > 0 )
            {
               SCIPhashmapFree(&conssmap);
            }
            if( SCIPgetNVars(scip) > 0 )
            {
               SCIPhashmapFree(&varmap);
            }
            SCIP_CALL_ABORT( SCIPtransformProb(paraInstanceScip));
            nCopies = 2;
            std::cout << "** ParaScipInstance is copied twice. **" << std::endl;
         }

         nVars = SCIPgetNVars(paraInstanceScip);
         SCIP_VAR **vars = SCIPgetVars(paraInstanceScip);

         /* make varName and objCoefs and ovnm */
         posVarNames = new int[nVars];
         objCoefs = new SCIP_Real[nVars];
         // mapToOriginalIndecies = new int[nVars];

         lVarNames = 0;
         for(int v = 0; v < nVars; ++v)
         {
            posVarNames[v] = lVarNames;
            objCoefs[v] = SCIPvarGetObj(vars[v]);
            assert(SCIPvarGetProbindex(vars[v])!=-1);
            assert(SCIPvarGetProbindex(vars[v]) == v);
            lVarNames += strlen(SCIPvarGetName(vars[v])) + 1;
         }
         varNames = new char[lVarNames];
         varLbs = new SCIP_Real[nVars];
         varUbs = new SCIP_Real[nVars];
         varTypes = new int[nVars];
         for(int v = 0; v < nVars; ++v )
         {
            SCIP_VAR *var = vars[v];
            strcpy (&varNames[posVarNames[v]], SCIPvarGetName(var) );
            varLbs[SCIPvarGetProbindex(var)] = SCIPvarGetLbLocal(var); //* we should use global?
            varUbs[SCIPvarGetProbindex(var)] = SCIPvarGetUbLocal(var); //* we should use global?
            varTypes[SCIPvarGetProbindex(var)] = SCIPvarGetType(var);
         }
      }
   }

   /** destractor */
   ~ScipParaInstanceMpi(
	        )
   {
      if( orgScip )
      {
         SCIPfree(&paraInstanceScip);
      }
   }

  /** create presolved problem instance that is solved by ParaSCIP form scip environment in this object */
  void copyScipEnvironment(
        SCIP **scip
        )
  {
     /** this routine for Pthread version. So, this should not be used **/
     abort();
  }

  SCIP *getScip(
        )
  {
     /** this routine for Pthread version. So, this should not be used **/
     abort();
  }

  void setFileName(const char *inFileName)
  {
     fileName = inFileName;
  }

   /** broadcasts instance to all solvers */
   int bcast(UG::ParaComm *comm, int rank, int method);



};

typedef ScipParaInstanceMpi *ScipParaInstanceMpiPtr;

}

#endif  // __SCIP_PARA_INSTANCE_MPI_H__

