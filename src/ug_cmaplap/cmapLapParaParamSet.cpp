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

/**@file    cmapLapParaParamSet.cpp
 * @brief   Parameter set for CMAP-LAP.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cfloat>
#include <climits>
#include "cmapLapParaParamSet.h"
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{


CMapLapParaParamSet::CMapLapParaParamSet(
      )
#if defined(_COMM_PTH) || defined(_COMM_CPP11)
      : UG::ParaParamSetTh(CMapLapParaParamsSize)
#else
      : UG::ParaParamSetMpi(CMapLapParaParamsSize)
#endif
{
   ///
   /// bool params
   ///
   paraParams[WarmStartOnlyPool] = new UG::ParaParamBool(
         "WarmStartOnlyPool",
         "# LoadCoordinator reads only Basis and Vector pool in warm start.",
         false);
   paraParams[NoWaitNotificationId] = new UG::ParaParamBool(
         "NoWaitNotificationId",
         "# Solvers do not wait NotificationId from LoadCoordinator (LoadCoordinator does not send recieve tag of NoWaitNotificationId)",
         false);
   paraParams[ShareIncumbentVector] = new UG::ParaParamBool(
         "ShareIncumbentVector",
         "# Each solver receives the incumbent solution [Default value: TRUE]",
         true);
   paraParams[LogShareDataPoolAll] = new UG::ParaParamBool(
         "LogShareDataPoolAll",
         "# Display norm of vectors in vector pool [Default value: FALSE]",
         false);
   paraParams[LogShareDataPoolStat] = new UG::ParaParamBool(
         "LogShareDataPoolStat",
         "# Display statistics information of vectors in vector pool [Default value: FALSE]",
         false);
   paraParams[CheckpointThreading] = new UG::ParaParamBool(
         "CheckpointThreading",
         "# Create a thread for writing checkpoint [Default value: FALSE]",
         true);
   paraParams[CheckpointReserving] = new UG::ParaParamBool(
         "CheckpointReserving",
         "# Copy reserved checkpoint objects. (This leads to increased memory usage.) [Default value: TRUE]",
         true);
   paraParams[LogSimilarityOfBasis] = new UG::ParaParamBool(
         "LogSimilarityOfBasis",
         "# LogSimilarityOfBasis is TRUE: output similarity of basis whose each solver [Default value: FALSE]",
         false);
   paraParams[DynamicDimensionOfSharedLattice] = new UG::ParaParamBool(
         "DynamicDimensionOfSharedLattice",
         "# Dimension of shared lattice will dynamically change in the session [Default value: FALSE]",
         false);
   paraParams[AutoAdjustmentNotificationInterval] = new UG::ParaParamBool(
         "AutoAdjustmentNotificationInterval",
         "# Adjust the solver's notification interval so that LCLowerIdleRatio <= idle ratio of LoadCoordinator <= LCUpperIdleRatio. [Default value: FALSE]",
         false);


   ///
   /// int params
   ///
   paraParams[MaxSizeOfMessageQueue] = new UG::ParaParamInt(
         "MaxSizeOfMessageQueue",
         "# Maximum size of the message queue of solver. Solver does not send any message other than paraSolution when the size of the message queue is greater than it. -1 means no limit of queue. [Default value: 100][-1, INT_MAX]",
         100,
         -1,
         INT_MAX);
   paraParams[BkzVerbose] = new UG::ParaParamInt(
         "BkzVerbose",
         "# Control log level of BKZ task: [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[EnumVerbose] = new UG::ParaParamInt(
         "EnumVerbose",
         "# Control log level of ENUM task: [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[SieveVerbose] = new UG::ParaParamInt(
         "SieveVerbose",
         "# Control log level of Sieve task: [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[NumRandomizedRows] = new UG::ParaParamInt(
         "NumRandomizedRows",
         "# Number of randomized rows from the bottom of basis. -1 means no randomization [Default value -1][-1, INT_MAX]",
         -1,
         -1,
         INT_MAX);
   paraParams[BkzStartBlockSize] = new UG::ParaParamInt(
         "BkzStartBlockSize",
         "# Bkz reduction with blocksize from BkzStartBlockSize to BkzEndBlockSize at intervals of BkzBlockSizeInterval [Default value: 30][2, INT_MAX]",
         30,
         2,
         INT_MAX);
   paraParams[BkzEndBlockSize] = new UG::ParaParamInt(
         "BkzEndBlockSize",
         "# Bkz reduction with blocksize from BkzStartBlockSize to BkzEndBlockSize at intervals of BkzBlockSizeInterval [Default value: 30][2, INT_MAX]",
         30,
         2,
         INT_MAX);
   paraParams[BkzBlockSizeInterval] = new UG::ParaParamInt(
         "BkzBlockSizeInterval",
         "# Bkz reduction with blocksize from BkzStartBlockSize to BkzEndBlockSize at intervals of BkzBlockSizeInterval [Default value: 30][1, INT_MAX]",
         5,
         1,
         INT_MAX);
   paraParams[BkzNumOfSendVectorsToPool] = new UG::ParaParamInt(
         "BkzNumOfSendVectorsToPool",
         "# Number of the vectors sended in Bkz to Load Coordnator [Default value: 1][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[BkzNumOfReceiveVectorsFromPool] = new UG::ParaParamInt(
         "BkzNumOfReceiveVectorsFromPool",
         "# Number of receiving vectors in Bkz from LoadCoordinator [Default value: 10][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[EnumNumOfSendVectorsToPool] = new UG::ParaParamInt(
         "EnumNumOfSendVectorsToPool",
         "# Number of sending vectors during Enumeration to LoadCoordinator [Default value: 1][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[EnumSamplingDepth] = new UG::ParaParamInt(
         "EnumSamplingDepth",
         "# Sampling depth during Enumeration [Default value: 2][0, INT_MAX]",
         2,
         0,
         INT_MAX);
   paraParams[EnumLowerNumberOfDividedSearchTree] = new UG::ParaParamInt(
         "EnumLowerNumberOfDividedSearchTree",
         "# Lower number of the divided enumeration trees. [Default value: 1000][1, INT_MAX]",
         1000,
         1,
         INT_MAX);
   paraParams[SubEnumProjectedDimension] = new UG::ParaParamInt(
         "SubEnumProjectedDimension",
         "# Projected dimension of SubEnum. 0 means run normal Enum [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[SubEnumNumOfSendVectorsToPool] = new UG::ParaParamInt(
         "SubEnumNumOfSendVectorsToPool",
         "# Number of sending vectors in SubEnum to LoadCoordinator [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[SieveMaxCollision] = new UG::ParaParamInt(
         "SieveMaxCollision",
         "# Sieve terminates the number of collision exceeds SieveMaxCollision. -1 means no limit [Default value: -1][-1, INT_MAX]",
         -1,
         -1,
         INT_MAX);
   paraParams[SieveMaxListSize] = new UG::ParaParamInt(
         "SieveMaxListSize",
         "# Maximum List L size. -1 means no limit [Default value: 10000000][-1, INT_MAX]",
         10000000,
         -1,
         INT_MAX);
   paraParams[SieveMaxStackSize] = new UG::ParaParamInt(
         "SieveMaxStackSize",
         "# Maximum Stack S size. -1 means no limit [Default value: 10000000][-1, INT_MAX]",
         10000000,
         -1,
         INT_MAX);
   paraParams[SieveStartNVectors] = new UG::ParaParamInt(
         "SieveStartNVectors",
         "# Sieve algorithm start the number of vectors in ShareDataPool is greater than SieveStartNVectors [Default value: 100][-1, INT_MAX]",
         100,
         -1,
         INT_MAX);
   paraParams[SieveNumThreads] = new UG::ParaParamInt(
         "SieveNumThreads",
         "# Maximum number of threads of a Sieve solver. 0 means the number of threads per rank in the machine [Default value: -1][-1, INT_MAX]",
         -1,
         -1,
         INT_MAX);
   paraParams[SieveNumOfSendVectorsToPool] = new UG::ParaParamInt(
         "SieveNumOfSendVectorsToPool",
         "# Number of sending vectors in Sieve to LoadCoordinator [Default value: 10][0, INT_MAX]",
         10,
         0,
         INT_MAX);
   paraParams[SieveNumOfReceiveVectorsFromPool] = new UG::ParaParamInt(
         "SieveNumOfReceiveVectorsFromPool",
         "# Number of receiving vectors in Sieve from LoadCoordinator [Default value: 100][0, INT_MAX]",
         100,
         0,
         INT_MAX);
   paraParams[BlockSizeForLocalSolver] = new UG::ParaParamInt(
         "BlockSizeForLocalSolver",
         "# blocksize of Bkz for Local Solver [Default value: 10][1, INT_MAX]",
         10,
         1,
         INT_MAX);
   paraParams[NumOfInitialBkzSolvers] = new UG::ParaParamInt(
         "NumOfInitialBkzSolvers",
         "# Number of Bkz solvers at the beginning of the run. -1 means it is number of solvers [Default value: -1][-1, INT_MAX]",
         -1,
         -1,
         INT_MAX);
   paraParams[NumOfInitialEnumSolvers] = new UG::ParaParamInt(
         "NumOfInitialEnumSolvers",
         "# Number of Enum solvers at the beginning of the run [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[NumOfInitialSieveSolvers] = new UG::ParaParamInt(
         "NumOfInitialSieveSolvers",
         "# Number of Sieve solvers at the beginning of the run [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[NumOfLocalSolversInLC] = new UG::ParaParamInt(
         "NumOfLocalSolversInLC",
         "# Number of local solvers in LoadCoordinator [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[ShareDataPoolSize] = new UG::ParaParamInt(
         "ShareDataPoolSize",
         "# Limit number of vectors stored in vector pool  [Default value: 100000][0, INT_MAX]",
         100000,
         0,
         INT_MAX);
   paraParams[InstancePoolSize] = new UG::ParaParamInt(
         "InstancePoolSize",
         "# Limit number of basis stored in basis pool  [Default value: 50][1, INT_MAX]",
         50,
         1,
         INT_MAX);
   paraParams[WriteSizeShareDataPool] = new UG::ParaParamInt(
         "WriteSizeShareDataPool",
         "# Maximum number of vector pool vectors that are written in every checkpoint  [Default value: 1000][1, INT_MAX]",
         1000,
         1,
         INT_MAX);
   paraParams[LocalTaskNVectorLowerBound] = new UG::ParaParamInt(
         "LocalTaskNVectorLowerBound",
         "# if number of vectors in vector pool is smaller than it, local task is not created [Default value: 10][2, INT_MAX]",
         10,
         2,
         INT_MAX);
   paraParams[LowerBoundOfInstancePoolShouldHold] = new UG::ParaParamInt(
         "LowerBoundOfInstancePoolShouldHold",
         "# Lower bound of number of basis that basis element pool tries to hold [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[ProjectedDimension] = new UG::ParaParamInt(
         "ProjectedDimension",
         "# projected dimension of this solver [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[DimensionOfSharedLattice] = new UG::ParaParamInt(
         "DimensionOfSharedLattice",
         "# dimension of always shared lattice [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[NumOfSamplingForSimilarityOfBasis] = new UG::ParaParamInt(
         "NumOfSamplingForSimilarityOfBasis",
         "# number of sampling for calculation of basis similarity [Default value: 100][2, INT_MAX]",
         100,
         2,
         INT_MAX);

   ///
   /// longint params
   ///

   ///
   /// real params
   ///
   paraParams[LowerBoundOfNorm] = new UG::ParaParamReal(
         "LowerBoundOfNorm",
         "# Solvers terminate when they finds the vector whose norm shorter than it. -1 means it is ignored [Default value: -1][0, DBL_MAX]",
         -1,
         0.0,
         DBL_MAX);
   paraParams[LowerBoundOfApproxFactor] = new UG::ParaParamReal(
         "LowerBoundOfApproxFactor",
         "# Solvers terminate when they finds the vector whose approx factor shorter than it. -1 means it is ignored [Default value: -1][0, DBL_MAX]",
         -1,
         0.0,
         DBL_MAX);
   paraParams[BkzTaskMaxTimeLimit] = new UG::ParaParamReal(
         "BkzTaskMaxTimeLimit",
         "# Max time limit of a Bkz para task. [Default value: DBL_MAX][0.0, DBL_MAX]",
         DBL_MAX,
         0.0,
         DBL_MAX);
   paraParams[EnumTaskMaxTimeLimit] = new UG::ParaParamReal(
         "EnumTaskMaxTimeLimit",
         "# Max time limit of a Enum para task. [Default value: DBL_MAX][0.0, DBL_MAX]",
         DBL_MAX,
         0.0,
         DBL_MAX);
   paraParams[SieveTaskMaxTimeLimit] = new UG::ParaParamReal(
         "SieveTaskMaxTimeLimit",
         "# Max time limit of a Sieve para task. [Default value: DBL_MAX][0.0, DBL_MAX]",
         DBL_MAX,
         0.0,
         DBL_MAX);
   paraParams[EnumTraversedNodesPerSeconds] = new UG::ParaParamReal(
         "EnumTraversedNodesPerSeconds",
         "# of nodes in an enumeration tree that can be traversed per second. [Default value: 3.3554432e7][1.0, DBL_MAX]",
         3.3554432e7,
         1.0,
         DBL_MAX);
   paraParams[EnumPruningParameter] = new UG::ParaParamReal(
         "EnumPruningParameter",
         "# pruning parameter of Enum Task. [Default value: 1.0][-1.0, 1.0]",
         1.0,
         -1.0,
         1.0);
   paraParams[RandomizeScale] = new UG::ParaParamReal(
         "RandomizeScale",
         "# When we generate randomize matrix, we multiplied RandomizeScale for each elements of matrix. [Default value: 1.2][1.0, DBL_MAX]",
         1.2,
         1.0,
         DBL_MAX);
   paraParams[IntervalTimeOfAssignmentTableOutput] = new UG::ParaParamReal(
         "IntervalTimeOfAssignmentTableOutput",
         "# Interval time of assignment table output [Default value: 5.0][0.0, DBL_MAX]",
         5.0,
         0.0,
         DBL_MAX);
   paraParams[IntervalTimeOfLogShareDataPool] = new UG::ParaParamReal(
         "IntervalTimeOfLogShareDataPool",
         "# Interval time to write log of vector pool [Default: 600.0][0.0, DBL_MAX]",
         600.0,
         0.0,
         DBL_MAX);
   paraParams[IReceiveInterval] = new UG::ParaParamReal(
         "IReceiveInterval",
         "# An active Solver calls iReceive() when time elapsed from its previous notification. [Default: 1.0][0.0, DBL_MAX]",
         1.0,
         0.0,
         DBL_MAX);
   paraParams[ShareVectorsInterval] = new UG::ParaParamReal(
         "ShareVectorsInterval",
         "# An active Solver sends vectors when time elapsed from its previous sending. [Default: 1.0][0.0, DBL_MAX]",
         1.0,
         0.0,
         DBL_MAX);
   paraParams[IntervalTimeOfLogSimilarityOfBasis] = new UG::ParaParamReal(
         "IntervalTimeOfLogSimilarityOfBasis",
         "# Interval time to output the similarity of basis whose each solver. [Default: 600.0][0.0, DBL_MAX]",
         600.0,
         0.0,
         DBL_MAX);
   paraParams[LCTermTimeUpdateNotificationInterval] = new UG::ParaParamReal(
         "LCTermTimeUpdateNotificationInterval",
         "# Interval time to update the solver's NotificationInterval. [Default: 600.0][0.0, DBL_MAX]",
         600.0,
         0.0,
         DBL_MAX);
   paraParams[LCUpperIdleRatio] = new UG::ParaParamReal(
         "LCUpperIdleRatio",
         "# Upperbound of ratio of LC's idle. [Default: 1.0][0.0, 1.0]",
         1.0,
         0.0,
         1.0);
   paraParams[LCLowerIdleRatio] = new UG::ParaParamReal(
         "LCLowerIdleRatio",
         "# Lowerbound of ratio of LC's idle. [Default: 0.3][0.0, 1.0]",
         0.3,
         0.0,
         1.0);



   ///
   /// char params
   ///

   ///
   /// string params
   ///
   paraParams[LaptoolsParamFilePath] = new UG::ParaParamString(
         "LaptoolsParamFilePath",
         "# Laptools parameter settings file name. Empty name use default settings. [Default: ]",
         "");

}

void
CMapLapParaParamSet::read(
      UG::ParaComm *comm,
      const char* filename
      )
{

   ParaParamSet::read(comm, filename);

   ///
   /// check parameter consistency
   ///
   if( (getIntParamValue(NumOfInitialBkzSolvers) < 0)
         + (getIntParamValue(NumOfInitialEnumSolvers) < 0)
         + (getIntParamValue(NumOfInitialSieveSolvers) < 0) >= 2 )
   {
      THROW_LOGICAL_ERROR7(
            "Must not be more than two negative values among NumOfInitialBkzSolvers, NumOfInitialEnumSolvers, and NumOfInitialSieveSolvers",
            "; NumOfInitialBkzSolvers = ",
            getIntParamValue(NumOfInitialBkzSolvers),
            "; NumOfInitialEnumSolvers = ",
            getIntParamValue(NumOfInitialEnumSolvers),
            "; NumOfInitialSieveSolvers = ",
            getIntParamValue(NumOfInitialSieveSolvers)
            );
   }

   if( getBoolParamValue(UG::Quiet) )
   {
      setBoolParamValue(UG::TagTrace, false);
      setBoolParamValue(UG::LogSolvingStatus, false);
      setBoolParamValue(UG::LogTasksTransfer, false);
   }

}

} // namespace ParaCMapLAP
