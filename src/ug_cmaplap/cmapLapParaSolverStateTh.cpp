/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*  CMAP-LAP --- Configurable Massively Parallel Solver for Lattice Problems */
/*                                                                           */
/*  Copyright Written by Nariaki Tateiwa <n-tateiwa@kyudai.jp>,              */
/*                       Yuji Shinano <shinano@zib.de>,                      */
/*            Copyright (C) 2021 by Zuse Institute Berlin,                   */
/*            licensed under LGPL version 3 or later.                        */
/*            Commercial licenses are available through <licenses@zib.de>    */
/*                                                                           */
/* This code is free software; you can redistribute it and/or                */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.     */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    cmapLapParaSolverStateTh.cpp
 * @brief   CMapLapParaSolverState extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaSolverStateTh.h"
#ifdef _COMM_MPI_WORLD
#include "cmapLapParaSolverLocalComm.h"
#endif
#include <assert.h>
#include <memory>
#include "ug/paraComm.h"
#include "ug/paraDef.h"
#include "ug/paraTagDef.h"
#include "cmapLapParaComm.h"
#include "cmapLapParaSolverState.h"
#include "cmapLapParaDef.h"


namespace ParaCMapLAP
{


CMapLapParaSolverStateTh*
CMapLapParaSolverStateTh::createDatatype(
      )
{
   switch( solverType )
   {
   case Bkz:
   {
      int *basis = 0;
      if( dimension * basisRows > 0 )
      {
         basis = new int [ dimension * basisRows ];
         for( int i = 0; i < dimension * basisRows; i++ )
         {
            basis[i] = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.basis[i];
         }
      }
      return new CMapLapParaSolverStateTh(
            notificationId,
            lcId,
            globalSubtreeIdInLc,
            detTime,
            threadId,
            dimension,
            basisRows,
            meanMessageQueueSize,
            maxMessageQueueSize,
            basis,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.currentBlockSize,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.tour,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.elapsedTime,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.shortestNorm,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.approxFactor,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.hermiteFactor,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.rootHermiteFactor,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCost,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCostGH,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.slopeGSA,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.topHalfSlopeGSA,
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.orthogonalFactor
            );
   }
   case Enum:
   {
      assert( basisRows > 0 );
      int* coeffs = new int [ basisRows ];
      for(int i = 0; i < basisRows; i++)
      {
         coeffs[i] = cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs[i];
      }
      return new CMapLapParaSolverStateTh(
            notificationId,
            lcId,
            globalSubtreeIdInLc,
            detTime,
            threadId,
            dimension,
            basisRows,
            meanMessageQueueSize,
            maxMessageQueueSize,
            cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.elapsedTime,
            coeffs,
            cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.searchNodeIndex,
            cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.enumCost,
            cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.shortestNorm,
            cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.approxFactor,
            cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.hermiteFactor,
            cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.rootHermiteFactor,
            cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.numSearchedNodes
            );
   }
   case Sieve:
   {
      return new CMapLapParaSolverStateTh(
            notificationId,
            lcId,
            globalSubtreeIdInLc,
            detTime,
            threadId,
            dimension,
            basisRows,
            meanMessageQueueSize,
            maxMessageQueueSize,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.elapsedTime,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.blockSize,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nLoop,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.listSize,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.stackSize,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.maxListSize,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nCollisions,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.shortestNorm,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.approxFactor,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.hermiteFactor,
            cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.rootHermiteFactor
            );
   }
   default:
      THROW_LOGICAL_ERROR2("CMapLapParaSolverStateTh::createDatatype: Invalid solver type = ", static_cast<int>(solverType));
   }
}

void
CMapLapParaSolverStateTh::send(
      UG::ParaComm *comm,
      int destination,
      int tag
      )
{

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   PARA_COMM_CALL(
      commTh->uTypeSend((void *)createDatatype(), UG::ParaSolverStateType, destination, tag)
   );
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
      localComm->uTypeSend((void *)createDatatype(), ParaSolverStateType, destination, tag)
   );
#endif
}

void
CMapLapParaSolverStateTh::receive(
      UG::ParaComm *comm,
      int source,
      int tag
      )
{
   std::unique_ptr<CMapLapParaSolverStateTh> received;

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   PARA_COMM_CALL(
      commTh->uTypeReceive((void **)&received, UG::ParaSolverStateType, source, tag)
   );
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
      localComm->uTypeReceive((void **)&received, ParaSolverStateType, source, tag)
   );
#endif

   racingStage          = received->racingStage;
   notificationId       = received->notificationId;
   lcId                 = received->lcId;
   globalSubtreeIdInLc  = received->globalSubtreeIdInLc;
   detTime              = received->detTime;
   threadId             = received->threadId;
   dimension            = received->dimension;
   basisRows            = received->basisRows;
   meanMessageQueueSize = received->meanMessageQueueSize;
   maxMessageQueueSize  = received->maxMessageQueueSize;

   switch( received->solverType )
   {
   case Bkz:
      solverType = Bkz;
      if( dimension * basisRows > 0 )
      {
         cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.basis = new int [ dimension * basisRows ];
         for( int i = 0; i < dimension * basisRows; i++ )
         {
            cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.basis[i] = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.basis[i];
         }
      }
      else
      {
         cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.basis = 0;  // empty basis
      }
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.currentBlockSize  = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.currentBlockSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.tour              = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.tour;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.elapsedTime       = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.elapsedTime;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.shortestNorm      = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.shortestNorm;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.approxFactor      = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.approxFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.hermiteFactor     = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.hermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.rootHermiteFactor = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.rootHermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCost          = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCost;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCostGH        = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCostGH;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.slopeGSA          = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.slopeGSA;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.topHalfSlopeGSA   = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.topHalfSlopeGSA;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.orthogonalFactor  = received->cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.orthogonalFactor;
      break;
   case Enum:
      solverType = Enum;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.elapsedTime = received->cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.elapsedTime;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs = new int [basisRows];
      for(int i = 0; i < dimension; i++)
      {
         cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs[i] = received->cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs[i];
      }
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.searchNodeIndex  = received->cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.searchNodeIndex;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.enumCost          = received->cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.enumCost;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.shortestNorm      = received->cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.shortestNorm;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.approxFactor      = received->cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.approxFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.hermiteFactor     = received->cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.hermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.rootHermiteFactor = received->cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.rootHermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.numSearchedNodes  = received->cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.numSearchedNodes;
      break;
   case Sieve:
      solverType = Sieve;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.elapsedTime       = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.elapsedTime;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.blockSize         = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.blockSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nLoop             = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nLoop;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.listSize          = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.listSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.stackSize         = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.stackSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.maxListSize       = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.maxListSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nCollisions       = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nCollisions;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.shortestNorm      = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.shortestNorm;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.approxFactor      = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.approxFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.hermiteFactor     = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.hermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.rootHermiteFactor = received->cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.rootHermiteFactor;
      break;
   default:
      THROW_LOGICAL_ERROR2("CMapLapParaSolverStateTh::receive: Invalid solver type = ", static_cast<int>(solverType));
   }

}

} // namespace ParaCMapLAP
