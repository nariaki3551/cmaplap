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

/**@file    paraSolverState.h
 * @brief   Base class for solver state.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SOLVER_STATE_H__
#define __CMAP_LAP_PARA_SOLVER_STATE_H__

#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cassert>
#include <memory>
#include "ug/paraSolverState.h"
#include "cmapLapParaLog.h"

namespace ParaCMapLAP
{

struct CMapLapParaSolverStateBkz {
   int*     basis;                  ///< 1-dimension array of lattice basis row major
   int      currentBlockSize;       ///< current Bkz block size
   int      tour;                   ///< number of Bkz loop
   double   elapsedTime;            ///< elapsed time
   double   shortestNorm;           ///< the shortest norm found
   double   approxFactor;           ///< approximate factor
   double   hermiteFactor;          ///< hermite factor
   double   rootHermiteFactor;      ///< (hermite factor)^(1/dim)
   double   enumCost;               ///< approximated nodes of enumeration tree with incumbent radius
   double   enumCostGH;             ///< approximated nodes of enumeration tree with GH radius
   double   slopeGSA;               ///< slope of GSA line [Pool] better a than b when a > b
   double   topHalfSlopeGSA;        ///< slope of top-half GSA line [Pool] better a than b when a > b
   double   orthogonalFactor;       ///< orthogonal factor [Pool] better a than b when a < b
   int      nNotChangeRows;         ///< number of not-changing rows
};

struct CMapLapParaSolverStateEnum {
   double   elapsedTime;            ///< elapsed time
   int      *coeffs;                ///< size is dimension, coefficient of a current search node
   int      searchNodeIndex;        ///< current index of searched node in enumeration tree
   double   enumCost;               ///< approximated nodes of enumeration tree with incumbent radius
   double   shortestNorm;           ///< the shortest norm found
   double   approxFactor;           ///< approximate factor
   double   hermiteFactor;          ///< hermite factor
   double   rootHermiteFactor;      ///< (hermite false)^(1/dim)
   long int numSearchedNodes;       ///< number of searched nodes in the enumeration tree [Pool] better a than b when b > a
};

struct CMapLapParaSolverStateSieve {
   double      elapsedTime;         ///< elapsed time
   int         blockSize;           ///< blocksize
   long int    nLoop;               ///< number of Sieve algorithm loop
   int         listSize;            ///< size of List L
   int         stackSize;           ///< size of Stack S
   int         maxListSize;         ///< maximum size of List L up to the point
   int         nCollisions;         ///< number of collision
   double      shortestNorm;        ///< the shortest norm found
   double      approxFactor;        ///< approximate factor
   double      hermiteFactor;       ///< hermite factor
   double      rootHermiteFactor;   ///< (hermite false)^(1/dim)
};

union CMapLapParaSolverStateData {
   CMapLapParaSolverStateBkz     cmapLapParaSolverStateBkz;
   CMapLapParaSolverStateEnum    cmapLapParaSolverStateEnum;
   CMapLapParaSolverStateSieve   cmapLapParaSolverStateSieve;
};

///
/// class CMapLapParaSolverState
/// (ParaSolver state object for notification message)
///
class CMapLapParaSolverState : public UG::ParaSolverState, public std::enable_shared_from_this<CMapLapParaSolverState>
{

protected:

   int threadId;              ///< threadId = -1 when solverType is Sieve
   SolverType solverType;     ///< solver type of  state
   int dimension;             ///< the dimension of basis
   int basisRows;             ///< the number of basis vectors of the lattice basis
   int meanMessageQueueSize;  ///< mean of the message queue size
   int maxMessageQueueSize;   ///< mean of the message queue size
   CMapLapParaSolverStateData cmapLapParaSolverStateData;

public:

   ///
   /// default constructor
   ///
   CMapLapParaSolverState(
          )
          : UG::ParaSolverState(),
            threadId(-1),
            solverType(Undefined),
            dimension(-1),
            basisRows(-1)
   {
   }

   ///
   /// copy constructor
   ///
   CMapLapParaSolverState(
         const CMapLapParaSolverState& cmapLapParaSolverState
          )
          : UG::ParaSolverState(cmapLapParaSolverState)
   {
      threadId = cmapLapParaSolverState.threadId;
      solverType = cmapLapParaSolverState.solverType;
      dimension = cmapLapParaSolverState.dimension;
      basisRows = cmapLapParaSolverState.basisRows;
      cmapLapParaSolverStateData = cmapLapParaSolverState.cmapLapParaSolverStateData;
   }

   ///
   /// Constructor of Bkz
   ///
   CMapLapParaSolverState(
         unsigned int  inNotificationId,        ///< id for this notification
         int           inLcId,                  ///< lc id of current ParaNode
         int           inGlobalSubtreeId,       ///< global subtree id of current ParaNode
         double        inDetTime,               ///< deterministic time, -1: should be non-deterministic
         int           inThreadId,              ///< thread id
         int           inDimension,             ///< the dimension of SVP instance
         int           inBasisRows,             ///< the number of basis vectors of the lattice basis
         int           inMeanMessageQueueSize,  ///< mean of the message queue size
         int           inMaxMessageQueueSize,   ///< max of the message queue size
         int*          inBasis,                 ///< 1-dimension array
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
         : UG::ParaSolverState(0, inNotificationId, inLcId, inGlobalSubtreeId, inDetTime),
           threadId(inThreadId),
           solverType(Bkz),
           dimension(inDimension),
           basisRows(inBasisRows),
           meanMessageQueueSize(inMeanMessageQueueSize),
           maxMessageQueueSize(inMaxMessageQueueSize)
   {
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.basis = inBasis;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.currentBlockSize = inCurrentBlockSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.tour = inTour;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.elapsedTime = inElapsedTime;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.shortestNorm = inShortestNorm;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.approxFactor = inApproxFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.hermiteFactor = inHermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.rootHermiteFactor = inRootHermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCost = inEnumCost;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCostGH = inEnumCostGH;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.slopeGSA = inSlopeGSA;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.topHalfSlopeGSA = inTopHalfSlopeGSA;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.orthogonalFactor = inOrthogonalFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.nNotChangeRows = -1;
   }

   ///
   /// Constructor of Enum
   ///
   CMapLapParaSolverState(
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
         int*          inCoeffs,                ///< size is dimension, coefficient of a current search node
         int           inSearchNodeIndex,       ///< current index of searched node in enumeration tree
         double        inEnumCost,              ///< log of approximted nodes of enumeration tree with incumbent radius
         double        inShortestNorm,          ///< the shortest norm found
         double        inApproxFactor,          ///< approximate factor
         double        inHermiteFactor,         ///< hermite factor
         double        inRootHermiteFactor,     ///< (hermite false)^(1/dim)
         long int      inNumSearchedNodes       ///< number of searched nodes in the enumeration tree
         )
         : UG::ParaSolverState(0, inNotificationId, inLcId, inGlobalSubtreeId, inDetTime),
           threadId(inThreadId),
           solverType(Enum),
           dimension(inDimension),
           basisRows(inBasisRows),
           meanMessageQueueSize(inMeanMessageQueueSize),
           maxMessageQueueSize(inMaxMessageQueueSize)
   {
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.elapsedTime = inElapsedTime;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs = inCoeffs;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.searchNodeIndex = inSearchNodeIndex;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.enumCost = inEnumCost;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.shortestNorm = inShortestNorm;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.approxFactor = inApproxFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.hermiteFactor = inHermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.rootHermiteFactor = inRootHermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.numSearchedNodes = inNumSearchedNodes;
   }

   /// constructor of Sieve
   CMapLapParaSolverState(
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
         : UG::ParaSolverState(0, inNotificationId, inLcId, inGlobalSubtreeId, inDetTime),
           threadId(inThreadId),
           solverType(Sieve),
           dimension(inDimension),
           basisRows(inBasisRows),
           meanMessageQueueSize(inMeanMessageQueueSize),
           maxMessageQueueSize(inMaxMessageQueueSize)
   {
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.elapsedTime = inElapsedTime;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.blockSize = inBlockSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nLoop = inNLoop;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.listSize = inListSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.stackSize = inStackSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.maxListSize = inMaxListSize;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nCollisions = inNCollisions;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.shortestNorm = inShortestNorm;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.approxFactor = inApproxFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.hermiteFactor = inHermiteFactor;
      cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.rootHermiteFactor = inRootHermiteFactor;
   }


   ///
   /// destractor
   ///
   virtual ~CMapLapParaSolverState(
         )
   {
      switch( solverType )
      {
      case Bkz:
         if( cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.basis )
         {
            delete[] cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.basis;
         }
         break;
      case Enum:
         if( cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs )
         {
            delete[] cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs;
         }
         break;
      case Sieve:
         break;
      default:
         std::cerr << "CMapLapParaSolverState::~CMapLapParaSolverState Invarid sovler type = " << static_cast<int>(solverType) << std::endl;
         abort();
      }
   }

   ///
   /// get thread ID
   /// @return thread ID
   ///
   virtual int getThreadId(
         )
   {
      return threadId;
   }

   ///
   /// set thread ID
   /// @return thread ID
   ///
   virtual void setThreadId(
         int inThreadId
         )
   {
      threadId = inThreadId;
   }

   ///
   /// get solver type
   /// @return solver type
   ///
   virtual SolverType getSolverType(
         )
   {
      return solverType;
   }

   ///
   /// getter of dimension
   /// @return dimension of basis
   ///
   virtual int getDimension(
         )
   {
      return dimension;
   }

   ///
   /// getter of basisRows
   /// @return number of rows of basis
   ///
   virtual int getBasisRows(
         )
   {
      return basisRows;
   }

   ///
   /// get elapsed time
   /// @return elapsed time
   ///
   virtual int getElapsedTime(
         )
   {
      switch( solverType )
      {
      case Bkz:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.elapsedTime;
      case Enum:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.elapsedTime;
      case Sieve:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.elapsedTime;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverState::getElapsedTime, invalid solver type = ", static_cast<int>(solverType));
      }
   }

   ///
   /// get basis
   /// @return current basis
   ///
   virtual int* getBasis(
         )
   {
      assert( solverType == Bkz );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.basis;
   }

   ///
   /// get current block size
   /// @return current block size
   ///
   virtual int getCurrentBlockSize(
         )
   {
      assert( solverType == Bkz || solverType == Sieve );
      if( solverType == Bkz )
      {
         return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.currentBlockSize;
      }
      if( solverType == Sieve )
      {
         return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.blockSize;
      }
      return 0;
   }

   ///
   /// get tour
   /// @return tour
   ///
   virtual int getCurrentTour(
         )
   {
      assert( solverType == Bkz );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.tour;
   }

   ///
   /// get shortestNorm
   /// @return shortest vector norm
   ///
   virtual double getShortestNorm(
         )
   {
      switch( solverType )
      {
      case Bkz:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.shortestNorm;
      case Enum:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.shortestNorm;
      case Sieve:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.shortestNorm;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverState::getShortestNorm, invalid solver type = ", static_cast<int>(solverType));
      }
   }

   ///
   /// get approxFactor
   /// @return approximation factor
   ///
   virtual double getApproxFactor(
         )
   {
      switch( solverType )
      {
      case Bkz:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.approxFactor;
      case Enum:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.approxFactor;
      case Sieve:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.approxFactor;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverState::getApproxFactor, invalid solver type = ", static_cast<int>(solverType));
      }
   }

   ///
   /// get hermiteFactor
   /// @return hermite factor
   ///
   virtual double getHermiteFactor(
         )
   {
      switch( solverType )
      {
      case Bkz:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.hermiteFactor;
      case Enum:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.hermiteFactor;
      case Sieve:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.hermiteFactor;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverState::getHermiteFactor, invalid solver type = ", static_cast<int>(solverType));
      }
   }

   ///
   /// get rootHermiteFactor
   /// @return hermite factor to 1/dimension power
   ///
   virtual double getRootHermiteFactor(
         )
   {
      switch( solverType )
      {
      case Bkz:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.rootHermiteFactor;
      case Enum:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.rootHermiteFactor;
      case Sieve:
         return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.rootHermiteFactor;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverState::getRootHermiteFactor, invalid solver type = ", static_cast<int>(solverType));
      }
   }

   ///
   /// get enumCost
   /// @return enumeration cost
   ///
   virtual double getEnumCost(
         )
   {
      assert( solverType == Bkz || solverType == Enum );
      if( solverType == Bkz )
         return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCost;
      else
         return cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.enumCost;
   }

   ///
   /// get slopeGSA
   /// @return the slope of GSA curve
   ///
   virtual double getSlopeGSA(
         )
   {
      assert( solverType == Bkz );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.slopeGSA;
   }

   ///
   /// get topHalfSlopeGSA
   /// @return the slope of GSA curve of top-half GSO vectors
   ///
   virtual double getTopHalfSlopeGSA(
         )
   {
      assert( solverType == Bkz );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.topHalfSlopeGSA;
   }

   ///
   /// get orthogonalFactor
   /// @return the orthogonal factor
   ///
   virtual double getOrthogonalFactor(
         )
   {
      assert( solverType == Bkz );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.orthogonalFactor;
   }

   ///
   /// get numSearchedNodes
   /// @return the number of nodes searched in Enum
   ///
   virtual long int getNumSearchedNodes(
         )
   {
      assert( solverType == Enum );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.numSearchedNodes;
   }

   ///
   /// get cloned coefficient
   /// @return *coeff
   ///
   virtual int * getClonedCoeffs(
         )
   {
      assert( solverType == Enum );
      int* coeffs = new int [dimension];
      for( int i = 0; i < dimension; i++)
      {
         coeffs[i] = cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs[i];
      }
      return coeffs;
   }

   ///
   /// get coefficient
   /// @return *coeff
   ///
   virtual int * getCoeffs(
         )
   {
      assert( solverType == Enum );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs;
   }

   ///
   /// get nLoop
   /// @return number of Sieve algorithm loop
   ///
   virtual unsigned int getNLoop(
         )
   {
      assert( solverType == Sieve );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nLoop;
   }

   ///
   /// getter of searchNodeIndex
   /// @return basis index of search node in enumeration
   ///
   virtual int getSearchNodeIndex(
         )
   {
      assert( solverType == Enum );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.searchNodeIndex;
   }

   ///
   /// get listSize
   /// @return list L size
   ///
   virtual int getListSize(
         )
   {
      assert( solverType == Sieve );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.listSize;
   }

   ///
   /// get stackSize
   /// @return stack S size
   ///
   virtual int getStackSize(
         )
   {
      assert( solverType == Sieve );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.stackSize;
   }

   ///
   /// get maxListSize
   /// @return maximum size of List L up to the point
   ///
   virtual int getMaxListSize(
         )
   {
      assert( solverType == Sieve );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.maxListSize;
   }

   ///
   /// get nCollisions
   /// @return number of collision is Sieve algorithm
   ///
   virtual int getNCollisions(
         )
   {
      assert( solverType == Sieve );
      return cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nCollisions;
   }

   ///
   /// getter of isRacingStage
   /// @return true if the Solver notified this message is in racing stage, false otherwise
   ///
   virtual bool isRacingStage(
         )
   {
      return (racingStage == 1);
   }

   ///
   /// getter of notification id
   /// @return notification id
   ///
   virtual unsigned int getNotificaionId(
         )
   {
      return notificationId;
   }

   ///
   /// getter of LoadCoordintor id
   /// @return LoadCoordinator id
   ///
   virtual int getLcId(
         )
   {
      return lcId;
   }

   ///
   /// getter of global subtree id
   /// @return global subtree id
   ///
   virtual int getGlobalSubtreeId(
         )
   {
      return globalSubtreeIdInLc;
   }

   ///
   /// getter of deterministic time
   /// @return deterministic time
   ///
   virtual double getDeterministicTime(
         )
   {
      return detTime;
   }

   ///
   /// get ObjectiveFunctionValue
   ///
   virtual double getObjectiveFunctionValue(
         )
   {
      switch( solverType )
      {
      case Bkz:
         return std::pow(cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.shortestNorm, 2);
         break;
      case Enum:
         return std::pow(cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.shortestNorm, 2);
         break;
      case Sieve:
         return std::pow(cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.shortestNorm, 2);
         break;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverPoolElement: Invalid solver type = ",static_cast<int>(solverType));
      }
   }

   ///
   /// set number of not-changing rows
   ///
   virtual void setNNotChangeRows(
         int nNotChangeRows
         )
   {
      assert( solverType == Bkz );
      cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.nNotChangeRows = nNotChangeRows;
   }

   ///
   /// stringfy CMapLapParaSolverState
   /// @return string to show inside of CMapLapParaSolverState
   ///
   virtual std::string toString(
         )
   {
      std::ostringstream s;
      s << "racingStage = " << racingStage << ", notificationId = " << notificationId << ": ";
      s << "[" << lcId << ":" << globalSubtreeIdInLc << "]";
      return s.str();
   }

   ///
   /// stringfy CMapLapParaSolverStateData
   /// @return string to show inside of CMapLapParaSolverStateData
   ///
   virtual std::string toStringLog(
         std::string delimiter=""
         )
   {
      switch( solverType )
      {
      case Bkz:
         return toStringLogBkz(delimiter);
      case Enum:
         return toStringLogEnum(delimiter);
      case Sieve:
         return toStringLogSieve(delimiter);
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverPoolElement: Invalid solver type = ",static_cast<int>(solverType));
      }
   }

   ///
   /// stringfy CMapLapParaSolverStateData
   /// @return string to show inside of CMapLapParaSolverStateBkz
   ///
   virtual const std::string toStringLogBkz(
         std::string delimiter=""
         )
   {
      std::string taskName       = "Bkz";
      double   elapsedTime       = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.elapsedTime;
      int      size              = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.currentBlockSize;
      long int iter              = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.tour;
      double   progress          = 0.0;
      double   leftTime          = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCostGH;
      double   logCost           = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.enumCost;
      double   shortestNorm      = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.shortestNorm;
      double   approxFactor      = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.approxFactor;
      double   hermiteFactor     = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.hermiteFactor;
      double   rootHermiteFactor = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.rootHermiteFactor;
      double   orthogonalFactor  = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.orthogonalFactor;
      int      nNotChangeRows    = cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.nNotChangeRows;
      std::ostringstream s_appendix;
      s_appendix << " -- rho " << std::fixed << std::setprecision(5)
                 << cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.slopeGSA
                 << ", half rho " << std::fixed << std::setprecision(5)
                 << cmapLapParaSolverStateData.cmapLapParaSolverStateBkz.topHalfSlopeGSA;
      s_appendix << ", not_change " << nNotChangeRows
                 << ", Mean " << meanMessageQueueSize;
                 // << ", Max "  << maxMessageQueueSize;

      return Logging::toStringLogBase(taskName, elapsedTime, size, iter, progress, leftTime, logCost,
            shortestNorm, approxFactor, hermiteFactor, rootHermiteFactor, orthogonalFactor,
            meanMessageQueueSize, maxMessageQueueSize,
            s_appendix.str(), delimiter);
   }

   ///
   /// stringfy CMapLapParaSolverStateData
   /// @return string to show inside of CMapLapParaSolverStateEnum
   ///
   virtual const std::string toStringLogEnum(
         std::string delimiter=""
         )
   {
      std::string taskName       = "Enum";
      double   elapsedTime       = cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.elapsedTime;
      int      size              = 0;
      long int iter              = cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.numSearchedNodes;
      double   logCost           = cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.enumCost;
      double   progress          = 0.0;
      double   leftTime          = 0.0;
      double   shortestNorm      = cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.shortestNorm;
      double   approxFactor      = cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.approxFactor;
      double   hermiteFactor     = cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.hermiteFactor;
      double   rootHermiteFactor = cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.rootHermiteFactor;
      double   orthogonalFactor  = -1;
      std::ostringstream s_appendix;
      s_appendix << " -- Mean " << meanMessageQueueSize
                 << ", Max "  << maxMessageQueueSize;

      return Logging::toStringLogBase(taskName, elapsedTime, size, iter, progress, leftTime, logCost,
            shortestNorm, approxFactor, hermiteFactor, rootHermiteFactor, orthogonalFactor,
            meanMessageQueueSize, maxMessageQueueSize,
            s_appendix.str(), delimiter);
   }

   ///
   /// stringfy CMapLapParaSolverStateData
   /// @return string to show inside of CMapLapParaSolverStateSieve
   ///
   virtual const std::string toStringLogSieve(
         std::string delimiter=""
         )
   {
      std::string taskName       = "Sieve";
      double   elapsedTime       = cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.elapsedTime;
      int      size              = cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.blockSize;
      long int iter              = cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nLoop;
      double   progress          = 0.0;
      double   leftTime          = 0.0;
      double   logCost           = 0.0;
      double   shortestNorm      = cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.shortestNorm;
      double   approxFactor      = cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.approxFactor;
      double   hermiteFactor     = cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.hermiteFactor;
      double   rootHermiteFactor = cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.rootHermiteFactor;
      double   orthogonalFactor = -1;
      std::ostringstream s_appendix;
      s_appendix << " -- |L| "
                 << cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.listSize
                 << " ,|S| "
                 << cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.stackSize;
      s_appendix << ", Mean " << meanMessageQueueSize
                 << ", Max "  << maxMessageQueueSize;

      return Logging::toStringLogBase(taskName, elapsedTime, size, iter, progress, leftTime, logCost,
            shortestNorm, approxFactor, hermiteFactor, rootHermiteFactor, orthogonalFactor,
            meanMessageQueueSize, maxMessageQueueSize,
            s_appendix.str(), delimiter);
   }

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_SOLVER_STATE_H__
