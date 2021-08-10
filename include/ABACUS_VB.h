#ifndef ABACUS_VB_H
#define ABACUS_VB_H

#include <math.h>
#include <numeric>
#include <queue>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <omp.h>

#include <Eigen/Eigen>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_vector.h>

// The Lieb-Liniger model
#include "SystemInfo.h"
#include "LiebLinState.h"
#include "MatrixElementsLiebLin.h"

// Scanning
#include "ScanInfo.h"
#include "ScanState.h"
#include "Scan.h"

// NRG routines
#include "InfoTargetState.h"
#include "QuenchInfo.h"
#include "TimeEvolutionInfo.h"
#include "NRGInfo.h"
#include "NRG.h"

const std::complex<double> II(0.0, 1.0);  //  Shorthand for i     
const static Eigen::IOFormat CSVReduced(Eigen::StreamPrecision, 
    Eigen::DontAlignCols, ",","\n");

#endif
