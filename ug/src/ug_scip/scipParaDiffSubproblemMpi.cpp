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

/**@file    scipParaDiffSubproblemMpi.cpp
 * @brief   ScipParaDiffSubproblem extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <mpi.h>
#include "scipParaCommMpi.h"
#include "scipParaDiffSubproblemMpi.h"

using namespace UG;
using namespace ParaSCIP;

/** create ScipDiffSubproblemDatatype1 */
/************************************************
 * Currently, Datatype1 is not necessary.       *
 * I create this code for the future extension. *
 ************************************************/
MPI_Datatype
ScipParaDiffSubproblemMpi::createDatatype1(
      )
{

   int nBlocks = 0;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

#ifndef UG_DEBUG_SOLUTION
   int blockLengths[8];
   MPI_Aint displacements[8];
   MPI_Datatype types[8];
#else
   int blockLengths[9];
   MPI_Aint displacements[9];
   MPI_Datatype types[9];
#endif

   MPI_CALL(
      MPI_Get_address( &localInfoIncluded, &startAddress )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = 0;
   types[nBlocks] = MPI_INT;
   nBlocks++;

   MPI_CALL(
      MPI_Get_address( &nBoundChanges, &address )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = address - startAddress;
   types[nBlocks] = MPI_INT;
   nBlocks++;

   if( branchLinearConss )
   {
      nBranchLinearConss = branchLinearConss->nLinearConss;
   }
   MPI_CALL(
      MPI_Get_address( &nBranchLinearConss, &address )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = address - startAddress;
   types[nBlocks] = MPI_INT;
   nBlocks++;

   if( branchSetppcConss )
   {
      nBranchSetppcConss = branchSetppcConss->nSetppcConss;
   }
   MPI_CALL(
      MPI_Get_address( &nBranchSetppcConss, &address )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = address - startAddress;
   types[nBlocks] = MPI_INT;
   nBlocks++;

   if( linearConss )
   {
      nLinearConss = linearConss->nLinearConss;
   }
   MPI_CALL(
      MPI_Get_address( &nLinearConss, &address )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = address - startAddress;
   types[nBlocks] = MPI_INT;
   nBlocks++;

   if( boundDisjunctions )
   {
      nBoundDisjunctions = boundDisjunctions->nBoundDisjunctions;
   }
   MPI_CALL(
      MPI_Get_address( &nBoundDisjunctions, &address )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = address - startAddress;
   types[nBlocks] = MPI_INT;
   nBlocks++;

   if( varBranchStats )
   {
      nVarBranchStats = varBranchStats->nVarBranchStats;
   }
   MPI_CALL(
      MPI_Get_address( &nVarBranchStats, &address )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = address - startAddress;
   types[nBlocks] = MPI_INT;
   nBlocks++;

   if( varValues )
   {
      nVarValueVars = varValues->nVarValueVars;
   }
   MPI_CALL(
      MPI_Get_address( &nVarValueVars, &address )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = address - startAddress;
   types[nBlocks] = MPI_INT;
   nBlocks++;

#ifdef UG_DEBUG_SOLUTION
   MPI_CALL(
      MPI_Get_address( &includeOptimalSol, &address )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = address - startAddress;
   types[nBlocks] = MPI_INT;
   nBlocks++;
#endif

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;
}

/** create ScipDiffSubproblemDatatype */
MPI_Datatype
ScipParaDiffSubproblemMpi::createDatatype2(
      bool memAllocNecessary
      )
{
   assert( nBoundChanges != 0 || nBranchLinearConss != 0 ||  nBranchSetppcConss != 0 || nLinearConss != 0 || nBoundDisjunctions != 0 || nVarBranchStats != 0 || nVarValueVars != 0 );

   int nBlocks = 0;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

   int blockLengths[42];           // reserve maximum number of elements
   MPI_Aint displacements[42];     // reserve maximum number of elements
   MPI_Datatype types[42];         // reserve maximum number of elements

   /* this duplicate send of &localInfoIncluded is a dummy to get startAddress */
   MPI_CALL(
      MPI_Get_address( &localInfoIncluded, &startAddress )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = 0;
   types[nBlocks] = MPI_INT;
   nBlocks++;

   if( nBoundChanges > 0 )
   {
      if( memAllocNecessary )
      {
         indicesAmongSolvers = new int[nBoundChanges];
         branchBounds = new SCIP_Real[nBoundChanges];
         boundTypes = new SCIP_BOUNDTYPE[nBoundChanges];
      }

      MPI_CALL(
         MPI_Get_address( indicesAmongSolvers, &address )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nBoundChanges;
      types[nBlocks] = MPI_INT;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( branchBounds, &address )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nBoundChanges;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundTypes, &address )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nBoundChanges;
      types[nBlocks] = MPI_INT;
      nBlocks++;
   }

   if( nBranchLinearConss > 0 )
   {
      if( memAllocNecessary )
      {
         branchLinearConss = new ScipParaDiffSubproblemBranchLinearCons();
         branchLinearConss->nLinearConss = nBranchLinearConss;
         branchLinearConss->linearLhss = new SCIP_Real[nBranchLinearConss];
         branchLinearConss->linearRhss = new SCIP_Real[nBranchLinearConss];
         branchLinearConss->nLinearCoefs = new int[nBranchLinearConss];
      }
      MPI_CALL(
         MPI_Get_address( branchLinearConss->linearLhss, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = branchLinearConss->nLinearConss;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( branchLinearConss->linearRhss, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = branchLinearConss->nLinearConss;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( branchLinearConss->nLinearCoefs, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = branchLinearConss->nLinearConss;
      types[nBlocks] = MPI_INT;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( &(branchLinearConss->lConsNames), &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = 1;
      types[nBlocks] = MPI_INT;
      nBlocks++;
   }

   if( nBranchSetppcConss > 0 )
   {
      if( memAllocNecessary )
      {
         branchSetppcConss = new ScipParaDiffSubproblemBranchSetppcCons();
         branchSetppcConss->nSetppcConss = nBranchSetppcConss;
         branchSetppcConss->nSetppcVars = new int[nBranchSetppcConss];
         branchSetppcConss->setppcTypes = new int[nBranchSetppcConss];
      }
      MPI_CALL(
         MPI_Get_address( branchSetppcConss->nSetppcVars, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = branchSetppcConss->nSetppcConss;
      types[nBlocks] = MPI_INT;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( branchSetppcConss->setppcTypes, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = branchSetppcConss->nSetppcConss;
      types[nBlocks] = MPI_INT;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( &(branchSetppcConss->lConsNames), &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = 1;
      types[nBlocks] = MPI_INT;
      nBlocks++;
   }

   if( nLinearConss > 0 )
   {
      if( memAllocNecessary )
      {
         linearConss = new ScipParaDiffSubproblemLinearCons();
         linearConss->nLinearConss = nLinearConss;
         linearConss->linearLhss = new SCIP_Real[nLinearConss];
         linearConss->linearRhss = new SCIP_Real[nLinearConss];
         linearConss->nLinearCoefs = new int[nLinearConss];
      }
      MPI_CALL(
         MPI_Get_address( linearConss->linearLhss, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = linearConss->nLinearConss;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( linearConss->linearRhss, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = linearConss->nLinearConss;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( linearConss->nLinearCoefs, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = linearConss->nLinearConss;
      types[nBlocks] = MPI_INT;
      nBlocks++;
   }

   if( nBoundDisjunctions > 0 )
   {
      if( memAllocNecessary )
      {
         boundDisjunctions = new ScipParaDiffSubproblemBoundDisjunctions();
         boundDisjunctions->nBoundDisjunctions = nBoundDisjunctions;
         boundDisjunctions->nVarsBoundDisjunction = new int[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionInitial = new SCIP_Bool[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionSeparate = new SCIP_Bool[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionEnforce = new SCIP_Bool[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionCheck = new SCIP_Bool[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionPropagate = new SCIP_Bool[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionLocal = new SCIP_Bool[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionModifiable = new SCIP_Bool[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionDynamic = new SCIP_Bool[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionRemovable = new SCIP_Bool[nBoundDisjunctions];
         boundDisjunctions->flagBoundDisjunctionStickingatnode = new SCIP_Bool[nBoundDisjunctions];
      }
      MPI_CALL(
         MPI_Get_address( &(boundDisjunctions->nTotalVarsBoundDisjunctions), &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = 1;
      types[nBlocks] = MPI_INT;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->nVarsBoundDisjunction, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionInitial, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionSeparate, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionEnforce, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionCheck, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionPropagate, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionLocal, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionModifiable, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionDynamic, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionRemovable, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( boundDisjunctions->flagBoundDisjunctionStickingatnode, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = boundDisjunctions->nBoundDisjunctions;
      types[nBlocks] = MPI_UNSIGNED;
      nBlocks++;

   }

   if( nVarBranchStats > 0 )
   {
      if( memAllocNecessary )
      {
         varBranchStats = new ScipParaDiffSubproblemVarBranchStats();
         varBranchStats->nVarBranchStats = nVarBranchStats;
         varBranchStats->idxBranchStatsVars = new int[nVarBranchStats];
         varBranchStats->downpscost = new SCIP_Real[nVarBranchStats];
         varBranchStats->uppscost = new SCIP_Real[nVarBranchStats];
         varBranchStats->downvsids = new SCIP_Real[nVarBranchStats];
         varBranchStats->upvsids = new SCIP_Real[nVarBranchStats];
         varBranchStats->downconflen = new SCIP_Real[nVarBranchStats];
         varBranchStats->upconflen = new SCIP_Real[nVarBranchStats];
         varBranchStats->downinfer = new SCIP_Real[nVarBranchStats];
         varBranchStats->upinfer = new SCIP_Real[nVarBranchStats];
         varBranchStats->downcutoff = new SCIP_Real[nVarBranchStats];
         varBranchStats->upcutoff = new SCIP_Real[nVarBranchStats];
      }

      MPI_CALL(
         MPI_Get_address( &(varBranchStats->offset), &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = 1;
      types[nBlocks] = MPI_INT;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->idxBranchStatsVars, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_INT;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->downpscost, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->uppscost, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->downvsids, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->upvsids, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->downconflen, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->upconflen, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->downinfer, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->upinfer, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->downcutoff, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varBranchStats->upcutoff, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varBranchStats->nVarBranchStats;
      types[nBlocks] = MPI_DOUBLE;
      nBlocks++;

   }

   if( nVarValueVars > 0 )
   {
      if( memAllocNecessary )
      {
         varValues = new ScipParaDiffSubproblemVarValues();
         varValues->nVarValueVars = nVarValueVars;
         varValues->idxVarValueVars = new int[nVarValueVars];
         varValues->nVarValueValues = new int[nVarValueVars];
      }

      MPI_CALL(
         MPI_Get_address( &(varValues->nVarValues), &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = 1;
      types[nBlocks] = MPI_INT;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varValues->idxVarValueVars, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varValues->nVarValueVars;
      types[nBlocks] = MPI_INT;
      nBlocks++;

      MPI_CALL(
         MPI_Get_address( varValues->nVarValueValues, &address )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = varValues->nVarValueVars;
      types[nBlocks] = MPI_INT;
      nBlocks++;
   }

   assert(nBlocks <= 41);

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
   );

   return datatype;
}

/** create ScipDiffSubproblemDatatype */
MPI_Datatype
ScipParaDiffSubproblemMpi::createDatatype3(
      bool memAllocNecessary
      )
{
   assert( nBranchLinearConss != 0 || nBranchSetppcConss != 0 || nLinearConss != 0 || nBoundDisjunctions != 0 || nVarValueVars != 0 );

   int nBlocks = 0;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

   int nVarValueBlock = 0;
   if( nVarValueVars > 0 )
   {
      for(int i = 0; i <  varValues->nVarValueVars; i++ )
      {
         if(  varValues->nVarValueValues[i] > 0 )
         {
            nVarValueBlock += 9;
         }
      }
   }

   int nTotalBlocks = 1 + ((nBranchLinearConss*2) + 1) + (nBranchSetppcConss + 1) + (nLinearConss*2) + (nBoundDisjunctions*3) + nVarValueBlock;
   int *blockLengths = new int[nTotalBlocks];
   MPI_Aint *displacements = new MPI_Aint[nTotalBlocks];
   MPI_Datatype *types = new MPI_Datatype[nTotalBlocks];

   /* this duplicate send of &localInfoIncluded is a dummy to get startAddress */
   MPI_CALL(
      MPI_Get_address( &localInfoIncluded, &startAddress )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = 0;
   types[nBlocks] = MPI_INT;
   nBlocks++;

   if( nBranchLinearConss > 0 )
   {
      if( memAllocNecessary )
      {
         assert(branchLinearConss);
         branchLinearConss->linearCoefs = new SCIP_Real*[nBranchLinearConss];
         branchLinearConss->idxLinearCoefsVars = new int*[nBranchLinearConss];
         branchLinearConss->consNames = new char[branchLinearConss->lConsNames];
      }

      for(int i = 0; i < nBranchLinearConss; i++ )
      {
         if( memAllocNecessary )
         {
            branchLinearConss->linearCoefs[i] = new SCIP_Real[branchLinearConss->nLinearCoefs[i]];
            branchLinearConss->idxLinearCoefsVars[i] = new int[branchLinearConss->nLinearCoefs[i]];
         }
         MPI_CALL(
            MPI_Get_address( branchLinearConss->linearCoefs[i], &address )
         );
         displacements[nBlocks] = address - startAddress;
         blockLengths[nBlocks] = branchLinearConss->nLinearCoefs[i];
         types[nBlocks] = MPI_DOUBLE;
         nBlocks++;

         MPI_CALL(
            MPI_Get_address( branchLinearConss->idxLinearCoefsVars[i], &address )
         );
         displacements[nBlocks] = address - startAddress;
         blockLengths[nBlocks] = branchLinearConss->nLinearCoefs[i];
         types[nBlocks] = MPI_INT;
         nBlocks++;
      }
      MPI_CALL(
         MPI_Get_address( branchLinearConss->consNames, &address )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = branchLinearConss->lConsNames;
      types[nBlocks] = MPI_CHAR;
      nBlocks++;
   }

   if( nBranchSetppcConss > 0 )
   {
      if( memAllocNecessary )
       {
          assert(branchSetppcConss);
          branchSetppcConss->idxSetppcVars = new int*[nBranchSetppcConss];
          branchSetppcConss->consNames = new char[branchSetppcConss->lConsNames];
       }

       for(int i = 0; i < nBranchSetppcConss; i++ )
       {
          if( memAllocNecessary )
          {
             branchSetppcConss->idxSetppcVars[i] = new int[branchSetppcConss->nSetppcVars[i]];
          }
          MPI_CALL(
             MPI_Get_address( branchSetppcConss->idxSetppcVars[i], &address )
          );
          displacements[nBlocks] = address - startAddress;
          blockLengths[nBlocks] = branchSetppcConss->nSetppcVars[i];
          types[nBlocks] = MPI_INT;
          nBlocks++;
       }
       MPI_CALL(
          MPI_Get_address( branchSetppcConss->consNames, &address )
       );
       displacements[nBlocks] = address - startAddress;
       blockLengths[nBlocks] = branchSetppcConss->lConsNames;
       types[nBlocks] = MPI_CHAR;
       nBlocks++;
   }

   if( nLinearConss > 0 )
   {
      if( memAllocNecessary )
      {
         assert(linearConss);
         linearConss->linearCoefs = new SCIP_Real*[nLinearConss];
         linearConss->idxLinearCoefsVars = new int*[nLinearConss];
      }

      for(int i = 0; i < nLinearConss; i++ )
      {
         if( memAllocNecessary )
         {
            linearConss->linearCoefs[i] = new SCIP_Real[linearConss->nLinearCoefs[i]];
            linearConss->idxLinearCoefsVars[i] = new int[linearConss->nLinearCoefs[i]];
         }
         MPI_CALL(
            MPI_Get_address( linearConss->linearCoefs[i], &address )
         );
         displacements[nBlocks] = address - startAddress;
         blockLengths[nBlocks] = linearConss->nLinearCoefs[i];
         types[nBlocks] = MPI_DOUBLE;
         nBlocks++;

         MPI_CALL(
            MPI_Get_address( linearConss->idxLinearCoefsVars[i], &address )
         );
         displacements[nBlocks] = address - startAddress;
         blockLengths[nBlocks] = linearConss->nLinearCoefs[i];
         types[nBlocks] = MPI_INT;
         nBlocks++;
      }
   }

   if( nBoundDisjunctions > 0 )
   {
      if( memAllocNecessary )
      {
         assert( boundDisjunctions );
         boundDisjunctions->idxBoundDisjunctionVars = new int*[nBoundDisjunctions];
         boundDisjunctions->boundTypesBoundDisjunction = new SCIP_BOUNDTYPE*[nBoundDisjunctions];
         boundDisjunctions->boundsBoundDisjunction = new SCIP_Real*[nBoundDisjunctions];
      }

      for( int i = 0; i < nBoundDisjunctions; i++ )
      {
         if( memAllocNecessary )
         {
            boundDisjunctions->idxBoundDisjunctionVars[i] = new int[boundDisjunctions->nVarsBoundDisjunction[i]];
            boundDisjunctions->boundTypesBoundDisjunction[i] = new SCIP_BOUNDTYPE[boundDisjunctions->nVarsBoundDisjunction[i]];
            boundDisjunctions->boundsBoundDisjunction[i] = new SCIP_Real[boundDisjunctions->nVarsBoundDisjunction[i]];
         }
         MPI_CALL(
            MPI_Get_address( boundDisjunctions->idxBoundDisjunctionVars[i], &address )
         );
         displacements[nBlocks] = address - startAddress;
         blockLengths[nBlocks] = boundDisjunctions->nVarsBoundDisjunction[i];
         types[nBlocks] = MPI_INT;
         nBlocks++;

         MPI_CALL(
            MPI_Get_address( boundDisjunctions->boundTypesBoundDisjunction[i], &address )
         );
         displacements[nBlocks] = address - startAddress;
         blockLengths[nBlocks] = boundDisjunctions->nVarsBoundDisjunction[i];
         types[nBlocks] = MPI_UNSIGNED;   // actual SCIP_BoundType is enum
         nBlocks++;

         MPI_CALL(
            MPI_Get_address( boundDisjunctions->boundsBoundDisjunction[i], &address )
         );
         displacements[nBlocks] = address - startAddress;
         blockLengths[nBlocks] = boundDisjunctions->nVarsBoundDisjunction[i];
         types[nBlocks] = MPI_DOUBLE;
         nBlocks++;

      }
   }

   if( nVarValueVars > 0 )
   {
      if( memAllocNecessary )
      {
         assert(varValues);
         varValues->varValue            = new SCIP_Real*[nVarValueVars];
         varValues->varValueDownvsids   = new SCIP_Real*[nVarValueVars];
         varValues->varVlaueUpvsids     = new SCIP_Real*[nVarValueVars];
         varValues->varValueDownconflen = new SCIP_Real*[nVarValueVars];
         varValues->varValueUpconflen   = new SCIP_Real*[nVarValueVars];
         varValues->varValueDowninfer   = new SCIP_Real*[nVarValueVars];
         varValues->varValueUpinfer     = new SCIP_Real*[nVarValueVars];
         varValues->varValueDowncutoff  = new SCIP_Real*[nVarValueVars];
         varValues->varValueUpcutoff    = new SCIP_Real*[nVarValueVars];
      }

      for(int i = 0; i < nVarValueVars; i++ )
      {
         if(varValues-> nVarValueValues[i] > 0 )
         {
            if( memAllocNecessary )
            {
               varValues->varValue[i]            = new SCIP_Real[varValues->nVarValueValues[i]];
               varValues->varValueDownvsids[i]   = new SCIP_Real[varValues->nVarValueValues[i]];
               varValues->varVlaueUpvsids[i]     = new SCIP_Real[varValues->nVarValueValues[i]];
               varValues->varValueDownconflen[i] = new SCIP_Real[varValues->nVarValueValues[i]];
               varValues->varValueUpconflen[i]   = new SCIP_Real[varValues->nVarValueValues[i]];
               varValues->varValueDowninfer[i]   = new SCIP_Real[varValues->nVarValueValues[i]];
               varValues->varValueUpinfer[i]     = new SCIP_Real[varValues->nVarValueValues[i]];
               varValues->varValueDowncutoff[i]  = new SCIP_Real[varValues->nVarValueValues[i]];
               varValues->varValueUpcutoff[i]    = new SCIP_Real[varValues->nVarValueValues[i]];
            }

            MPI_CALL(
               MPI_Get_address( varValues->varValue[i], &address )
            );
            displacements[nBlocks] = address - startAddress;
            blockLengths[nBlocks] = varValues->nVarValueValues[i];
            types[nBlocks] = MPI_DOUBLE;
            nBlocks++;

            MPI_CALL(
               MPI_Get_address( varValues->varValueDownvsids[i], &address )
            );
            displacements[nBlocks] = address - startAddress;
            blockLengths[nBlocks] = varValues->nVarValueValues[i];
            types[nBlocks] = MPI_DOUBLE;
            nBlocks++;

            MPI_CALL(
               MPI_Get_address( varValues->varVlaueUpvsids[i], &address )
            );
            displacements[nBlocks] = address - startAddress;
            blockLengths[nBlocks] = varValues->nVarValueValues[i];
            types[nBlocks] = MPI_DOUBLE;
            nBlocks++;

            MPI_CALL(
               MPI_Get_address( varValues->varValueDownconflen[i], &address )
            );
            displacements[nBlocks] = address - startAddress;
            blockLengths[nBlocks] = varValues->nVarValueValues[i];
            types[nBlocks] = MPI_DOUBLE;
            nBlocks++;

            MPI_CALL(
               MPI_Get_address( varValues->varValueUpconflen[i], &address )
            );
            displacements[nBlocks] = address - startAddress;
            blockLengths[nBlocks] = varValues->nVarValueValues[i];
            types[nBlocks] = MPI_DOUBLE;
            nBlocks++;

            MPI_CALL(
               MPI_Get_address( varValues->varValueDowninfer[i], &address )
            );
            displacements[nBlocks] = address - startAddress;
            blockLengths[nBlocks] = varValues->nVarValueValues[i];
            types[nBlocks] = MPI_DOUBLE;
            nBlocks++;

            MPI_CALL(
               MPI_Get_address( varValues->varValueUpinfer[i], &address )
            );
            displacements[nBlocks] = address - startAddress;
            blockLengths[nBlocks] = varValues->nVarValueValues[i];
            types[nBlocks] = MPI_DOUBLE;
            nBlocks++;

            MPI_CALL(
               MPI_Get_address( varValues->varValueDowncutoff[i], &address )
            );
            displacements[nBlocks] = address - startAddress;
            blockLengths[nBlocks] = varValues->nVarValueValues[i];
            types[nBlocks] = MPI_DOUBLE;
            nBlocks++;

            MPI_CALL(
               MPI_Get_address( varValues->varValueUpcutoff[i], &address )
            );
            displacements[nBlocks] = address - startAddress;
            blockLengths[nBlocks] = varValues->nVarValueValues[i];
            types[nBlocks] = MPI_DOUBLE;
            nBlocks++;

         }
      }
   }


   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
   );

   delete [] blockLengths;
   delete [] displacements;
   delete [] types;

   return datatype;
}

int
ScipParaDiffSubproblemMpi::bcast(ParaComm *comm, int root)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI_Datatype datatype1;
   datatype1 = createDatatype1();
   MPI_CALL(
      MPI_Type_commit( &datatype1 )
   );
   PARA_COMM_CALL(
      commMpi->ubcast(&localInfoIncluded, 1, datatype1, root)
   );
   MPI_CALL(
      MPI_Type_free( &datatype1 )
   );

   if( nBoundChanges > 0 || nBranchLinearConss > 0 || nBranchSetppcConss > 0 || nLinearConss > 0 || nBoundDisjunctions > 0 || nVarBranchStats > 0 || nVarValueVars > 0 ){
      MPI_Datatype datatype2;
      if( comm->getRank() == root )
      {
         datatype2 = createDatatype2(false);
      }
      else
      {
         datatype2 = createDatatype2(true);
      }
      MPI_CALL(
         MPI_Type_commit( &datatype2 )
      );
      PARA_COMM_CALL(
            commMpi->ubcast(&localInfoIncluded, 1, datatype2, root)
      );

      MPI_CALL(
         MPI_Type_free( &datatype2 )
      );

      if( nBranchLinearConss > 0 || nBranchSetppcConss > 0 || nLinearConss > 0 || nBoundDisjunctions > 0 || nVarValueVars > 0 )
      {
         MPI_Datatype datatype3;
         if( comm->getRank() == root )
         {
            datatype3 = createDatatype3(false);
         }
         else
         {
            datatype3 = createDatatype3(true);
         }
         MPI_CALL(
            MPI_Type_commit( &datatype3 )
         );

         PARA_COMM_CALL(
               commMpi->ubcast(&localInfoIncluded, 1, datatype3, root)
         );
         MPI_CALL(
            MPI_Type_free( &datatype3 )
         );
      }
   }
   return 0;
}

int
ScipParaDiffSubproblemMpi::send(ParaComm *comm, int dest)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI_Datatype datatype1;
   datatype1 = createDatatype1();
   MPI_CALL(
      MPI_Type_commit( &datatype1 )
   );
   PARA_COMM_CALL(
      commMpi->usend(&localInfoIncluded, 1, datatype1, dest, TagDiffSubproblem)
   );
   MPI_CALL(
      MPI_Type_free( &datatype1 )
   );

   if( nBoundChanges > 0 || nBranchLinearConss > 0 || nBranchSetppcConss > 0 || nLinearConss > 0 || nBoundDisjunctions > 0 || nVarBranchStats > 0 || nVarValueVars > 0 )
   {
      MPI_Datatype datatype2;
      datatype2 = createDatatype2(false);
      MPI_CALL(
         MPI_Type_commit( &datatype2 )
      );
      PARA_COMM_CALL(
            commMpi->usend(&localInfoIncluded, 1, datatype2, dest, TagDiffSubproblem1)
      );
      MPI_CALL(
         MPI_Type_free( &datatype2 )
      );
      if( nBranchLinearConss > 0 || nBranchSetppcConss > 0 || nLinearConss > 0 || nBoundDisjunctions > 0  || nVarValueVars > 0 )
      {
         MPI_Datatype datatype3;
         datatype3 = createDatatype3(false);
         MPI_CALL(
            MPI_Type_commit( &datatype3 )
         );
         PARA_COMM_CALL(
               commMpi->usend(&localInfoIncluded, 1, datatype3, dest, TagDiffSubproblem2)
         );
         MPI_CALL(
            MPI_Type_free( &datatype3 )
         );
      }
   }
   return 0;
}

int
ScipParaDiffSubproblemMpi::receive(ParaComm *comm, int source)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI_Datatype datatype1;
   datatype1 = createDatatype1();
   MPI_CALL(
      MPI_Type_commit( &datatype1 )
   );
   PARA_COMM_CALL(
      commMpi->ureceive(&localInfoIncluded, 1, datatype1, source, TagDiffSubproblem)
   );
   MPI_CALL(
      MPI_Type_free( &datatype1 )
   );

   if( nBoundChanges > 0 || nBranchLinearConss > 0 || nBranchSetppcConss > 0 || nLinearConss > 0 || nBoundDisjunctions > 0 || nVarBranchStats > 0 || nVarValueVars > 0 )
   {
      MPI_Datatype datatype2;
      datatype2 = createDatatype2(true);
      MPI_CALL(
         MPI_Type_commit( &datatype2 )
      );
      PARA_COMM_CALL(
            commMpi->ureceive(&localInfoIncluded, 1, datatype2, source, TagDiffSubproblem1)
      );
      MPI_CALL(
         MPI_Type_free( &datatype2 )
      );
      if( nBranchLinearConss > 0 || nBranchSetppcConss > 0 || nLinearConss > 0 || nBoundDisjunctions > 0 || nVarValueVars > 0 )
      {
         MPI_Datatype datatype3;
         datatype3 = createDatatype3(true);
         MPI_CALL(
            MPI_Type_commit( &datatype3 )
         );
         PARA_COMM_CALL(
               commMpi->ureceive(&localInfoIncluded, 1, datatype3, source, TagDiffSubproblem2)
         );
         MPI_CALL(
            MPI_Type_free( &datatype3 )
         );
      }
   }
   return 0;
}

/** create clone of this object */
ScipParaDiffSubproblemMpi *
ScipParaDiffSubproblemMpi::clone(ParaComm *comm)
{
   return( new ScipParaDiffSubproblemMpi(this) );

}
