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

/**@file    cmapLapParaSolverStateMpi.h
 * @brief   CMapLapParaSolverState extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SOLVER_STATE_MPI_H__
#define __CMAP_LAP_PARA_SOLVER_STATE_MPI_H__


#include <memory>
#include <mpi.h>
#include "cmapLapParaSolverState.h"
#include "cmapLapParaDef.h"
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{

///
/// class CMapLapParaSolverStateMpi
/// (ParaSolver state object for notification message communicated by MPI)
///
class CMapLapParaSolverStateMpi : public CMapLapParaSolverState
{

   int iSolverType;          ///< int solver type: 0:Bkz, 1:Enum, 2:Sieve, -1: not set
   int typeAndDim[3];

   ///
   /// create Bkz CMapLapParaSolverStateMpi datatype
   /// @return MPI datatype for CMapLapParaSolverStateMpi
   ///
   virtual MPI_Datatype createDatatypeBkz(
         );

   ///
   /// create Enum CMapLapParaSolverStateMpi datatype
   /// @return MPI datatype for CMapLapParaSolverStateMpi
   ///
   virtual MPI_Datatype createDatatypeEnum(
         );

   ///
   /// create Sieve CMapLapParaSolverStateMpi datatype
   /// @return MPI datatype for CMapLapParaSolverStateMpi
   ///
   virtual MPI_Datatype createDatatypeSieve(
         );

public:

   ///
   /// default constructor
   ///
   CMapLapParaSolverStateMpi(
         )
         : iSolverType(static_cast<int>(Undefined))
   {
   }

   ///
   /// constructor Bkz
   ///
   CMapLapParaSolverStateMpi(
         unsigned int  inNotificationId,        ///< id for this notification
         int           inLcId,                  ///< lc id of current ParaNode
         int           inGlobalSubtreeId,       ///< global subtree id of current ParaNode
         double        inDetTime,               ///< deterministic time, -1: should be non-deterministic
         int           inThreadId,              ///< thread id
         int           inDimension,             ///< dimension of SVP instance
         int           inBasisRows,             ///< the number of basis vectors of the lattice basis
         int           inMeanMessageQueueSize,  ///< mean of the message queue size
         int           inMaxMessageQueueSize,   ///< max of the message queue size
         int*          inBasis,                 ///< 1-dimension array of lattice basis row major
         int           inCurrentBlockSize,      ///< current Bkz block size
         int           inTour,                  ///< number of Bkz loop
         double        inElapsedTime,           ///< elapsed time
         double        inShortestNorm,          ///< the shortest norm found
         double        inApproxFactor,          ///< approximated factor
         double        inHermiteFactor,         ///< hermite factor
         double        inRootHermiteFactor,     ///< (hermite factor)^(1/dim)
         double        inEnumCost,              ///< log of approximted nodes of enumeration tree with incumbent radius
         double        inEnumCostGH,            ///< log of approximted nodes of enumeration tree with GH radius
         double        inSlopeGSA,              ///< slope of GSA line
         double        inTopHalfSlopeGSA,       ///< slope of top-half GSA line
         double        inOrthogonalFactor       ///< orthogonal factor
         )
         : CMapLapParaSolverState(inNotificationId,
                                 inLcId,
                                 inGlobalSubtreeId,
                                 inDetTime,
                                 inThreadId,
                                 inDimension,
                                 inBasisRows,
                                 inMeanMessageQueueSize,
                                 inMaxMessageQueueSize,
                                 inBasis,
                                 inCurrentBlockSize,
                                 inTour,
                                 inElapsedTime,
                                 inShortestNorm,
                                 inApproxFactor,
                                 inHermiteFactor,
                                 inRootHermiteFactor,
                                 inEnumCost,
                                 inEnumCostGH,
                                 inSlopeGSA,
                                 inTopHalfSlopeGSA,
                                 inOrthogonalFactor
                                 ),
            iSolverType(static_cast<int>(Bkz))
   {
   }

   ///
   /// constructor Enum
   ///
   CMapLapParaSolverStateMpi(
         unsigned int  inNotificationId,        ///< id for this notification
         int           inLcId,                  ///< lc id of current ParaNode
         int           inGlobalSubtreeId,       ///< global subtree id of current ParaNode
         double        inDetTime,               ///< deterministic time, -1: should be non-deterministic
         int           inThreadId,              ///< thread id
         int           inDimension,             ///< dimension of SVP instance
         int           inBasisRows,             ///< the number of basis vectors of the lattice basis
         int           inMeanMessageQueueSize,  ///< mean of the message queue size
         int           inMaxMessageQueueSize,   ///< max of the message queue size
         double        inElapsedTime,           ///< elapsed time
         int*          inCoeff,                 ///< size is dimension, coefficient of a current search node
         int           inDepth,                 ///< current depth of searched node in enumeration tree
         double        inEnumCost,              ///< log of approximted nodes of enumeration tree with incumbent radius
         double        inShortestNorm,          ///< the shortest norm found
         double        inApproxFactor,          ///< approximate factor
         double        inHermiteFactor,         ///< hermite factor
         double        inRootHermiteFactor,     ///< (hermite false)^(1/dim)
         long int      inNumSearchedNodes       ///< number of searched nodes in the enumeration tree
         )
         : CMapLapParaSolverState(inNotificationId,
                                 inLcId,
                                 inGlobalSubtreeId,
                                 inDetTime,
                                 inThreadId,
                                 inDimension,
                                 inBasisRows,
                                 inMeanMessageQueueSize,
                                 inMaxMessageQueueSize,
                                 inElapsedTime,
                                 inCoeff,
                                 inDepth,
                                 inEnumCost,
                                 inShortestNorm,
                                 inApproxFactor,
                                 inHermiteFactor,
                                 inRootHermiteFactor,
                                 inNumSearchedNodes
                                 ),
            iSolverType(static_cast<int>(Enum))
   {
   }

   ///
   /// constructor of Sieve
   ///
   CMapLapParaSolverStateMpi(
         unsigned int  inNotificationId,        ///< id for this notification
         int           inLcId,                  ///< lc id of current ParaNode
         int           inGlobalSubtreeId,       ///< global subtree id of current ParaNode
         double        inDetTime,               ///< deterministic time, -1: should be non-deterministic
         int           inThreadId,              ///< thread id
         int           inDimension,             ///< dimension of SVP instance
         int           inBasisRows,             ///< the number of basis vectors of the lattice basis
         int           inMeanMessageQueueSize,  ///< mean of the message queue size
         int           inMaxMessageQueueSize,   ///< max of the message queue size
         double        inElapsedTime,           ///< elapsed time
         int           inBlockSize,             ///< block size
         long int      inNLoop,                 ///< number of Sieve algorithm loop
         int           inListSize,              ///< size of List L
         int           inStackSize,             ///< size of Stack S
         int           inMaxListSize,           ///< maximum size of List L up to the point
         int           inNCollisions,           ///< number of collision
         double        inShortestNorm,          ///< the shortest norm found
         double        inApproxFactor,          ///< approximated factor
         double        inHermiteFactor,         ///< hermite factor
         double        inRootHermiteFactor      ///< (hermite factor)^(1/dim)
         )
         : CMapLapParaSolverState(inNotificationId,
                                 inLcId,
                                 inGlobalSubtreeId,
                                 inDetTime,
                                 inThreadId,
                                 inDimension,
                                 inBasisRows,
                                 inMeanMessageQueueSize,
                                 inMaxMessageQueueSize,
                                 inElapsedTime,
                                 inBlockSize,
                                 inNLoop,
                                 inListSize,
                                 inStackSize,
                                 inMaxListSize,
                                 inNCollisions,
                                 inShortestNorm,
                                 inApproxFactor,
                                 inHermiteFactor,
                                 inRootHermiteFactor
                                 ),
            iSolverType(static_cast<int>(Sieve))
   {
   }

   ///
   /// destractor
   ///
   virtual ~CMapLapParaSolverStateMpi(
         )
   {
   }

   ///
   /// send this object
   /// @return always 0 (for future extensions)
   ///
   virtual void send(
         UG::ParaComm *comm,     ///< communicator used
         int destination,        ///< destination rank
         int tag                 ///< UG::TagSolverState
         );

   ///
   /// receive this object
   /// @return always 0 (for future extensions)
   ///
   virtual void receive(
         UG::ParaComm *comm,     ///< communicator used
         int source,             ///< source rank
         int tag                 ///< UG::TagSolverState
         );

};

#define DEF_PARA_SOLVER_STATE( para_state, state ) CMapLapParaSolverStateMpi *para_state = dynamic_cast< CMapLapParaSolverStateMpi* >(state)

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_SOLVER_STATE_MPI_H__
