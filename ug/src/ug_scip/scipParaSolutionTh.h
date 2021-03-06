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

/**@file    scipParaSolutionTh.h
 * @brief   ScipParaSolution extension for threads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_SOLUTION_TH_H__
#define __SCIP_PARA_SOLUTION_TH_H__

#include "ug/paraTagDefTh.h"
#ifdef _COMM_PTH
#include "ug/paraCommPth.h"
#endif
#ifdef _COMM_CPP11
#include "ug/paraCommCPP11.h"
#endif
#include "scipParaSolution.h"
#include "scipParaSolver.h"
#include "scip/scip.h"

namespace ParaSCIP
{

/** ScipSolution class */
class ScipParaSolutionTh : public ScipParaSolution
{
   /** create scipSolutionDatatype */
   ScipParaSolutionTh *createDatatype(UG::ParaComm *comm);

public:

   /** default constructor */
   ScipParaSolutionTh(
	      )
   {
   }

   /** constructor */
   ScipParaSolutionTh(
         ScipParaSolver *solver,
         SCIP_Real      objval,
         int            inNvars,
         SCIP_VAR **    vars,
         SCIP_Real *    vals
         )
        : ScipParaSolution(solver, objval, inNvars, vars, vals){}

   /** constructor */
   ScipParaSolutionTh(
         double inObjectiveFunctionValue,
         int inNVars,                       /**< number of variables */
         int *inIndicesAmongSolvers,        /**< array of variable indices ( probindex )  */
         SCIP_Real *inValues                /**< array of bounds which the branchings     */
         ): ScipParaSolution(inObjectiveFunctionValue, inNVars, inIndicesAmongSolvers, inValues) {}

   /** destructor */
   ~ScipParaSolutionTh(
         )
   {
   }

   /** create clone of this object */
   ScipParaSolutionTh *clone(UG::ParaComm *comm);

   /** broadcast solution data to from the root rank */
   void bcast(UG::ParaComm *comm, int root);

   /** send solution data to the rank */
   void send(UG::ParaComm *comm, int destination);

   /** receive solution data from the source rank */
   void receive(UG::ParaComm *comm, int source);

};

typedef ScipParaSolutionTh *ScipParaSolutionThPtr;

}

#endif // __SCIP_PARA_SOLUTION_TH_H__

