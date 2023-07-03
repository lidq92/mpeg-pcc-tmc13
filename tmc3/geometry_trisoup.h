/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2018, ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of the ISO/IEC nor the names of its contributors
 *   may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <cstdint>
#include <vector>
#include <cstring>

#include "PCCPointSet.h"
#include "geometry_octree.h"
#include "ringbuf.h"

namespace pcc {
//============================================================================

  struct CentroidDrift {
    int driftQ;
    int lowBound;
    int highBound;
    int ctxMinMax;
    int lowBoundSurface;
    int  highBoundSurface;

  };


  // Representation for a vertex in preparation for sorting.
  struct Vertex {
    Vec3<int32_t> pos;  // position of vertex
    int32_t theta;      // angle of vertex when projected along dominant axis
    int32_t tiebreaker;  // coordinate of vertex along dominant axis
    bool operator()(Vertex v1, Vertex v2)
    {
      if (v1.theta > v2.theta)
        return true;  // sort in decreasing order of theta
      if (v1.theta == v2.theta && v1.tiebreaker < v2.tiebreaker)
        return true;
      return false;
    }

    bool operator==(Vertex v1)
    {
      return this->pos == v1.pos;
    }

    bool operator<(Vertex v1)
    {
      return this->pos < v1.pos;
    }

  } ;


 //============================================================================
void determineTrisoupVertices(
  const ringbuf<PCCOctree3Node>& leaves,
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,  
  const PCCPointSet3& pointCloud,
  const int defaultBlockWidth,
  const int bitDropped,
  int distanceSearchEncoder);

void determineTrisoupNeighbours(
  const ringbuf<PCCOctree3Node>& leaves, 
  std::vector<uint16_t>& neighbNodes, 
  std::vector<int>& indexBefore,
  std::vector<std::vector<int>>& perpVertexStart,
  const int defaultBlockWidth);

void encodeTrisoupVertices(  
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<int>& indexBefore,
  std::vector<std::vector<int>>& perpVertexStart,
  int bitDropped,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  pcc::EntropyEncoder* arithmeticEncoder);

void decodeTrisoupVertices(  
  std::vector<bool>& segind,
  std::vector<uint8_t>& vertices,
  std::vector<uint16_t>& neighbNodes,
  std::vector<int>& indexBefore,
  std::vector<std::vector<int>>& perpVertexStart,
  int bitDropped,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  pcc::EntropyDecoder& arithmeticDecoder);


void encodeTrisoupCentroidResidue(
  std::vector<CentroidDrift>& drifts,
  pcc::EntropyEncoder* arithmeticEncoder);

void decodeTrisoupCommon(
  const ringbuf<PCCOctree3Node>& leaves,
  const std::vector<bool>& segind,
  const std::vector<uint8_t>& vertices,
  std::vector<CentroidDrift>& drifts,
  PCCPointSet3& pointCloud,
  PCCPointSet3& recPointCloud,
  int defaultBlockWidth,
  int poistionClipValue,
  uint32_t samplingValue,
  const int bitDropped,
  const bool isCentroidDriftActivated,
  bool isDecoder,
  bool haloFlag,
  bool adaptiveHaloFlag,
  bool fineRayflag,
  pcc::EntropyDecoder* arithmeticDecoder);

int findDominantAxis(
  std::vector<Vertex>& leafVertices,
  uint32_t blockWidth, 
  Vec3<int32_t> blockCentroid ); 

void rayTracingAlongdirection(
  std::vector<Vec3<int32_t>>& refinedVerticesBlock,
  int direction,
  uint32_t samplingValue,
  Vec3<int32_t> posNode,
  int minRange[3],
  int maxRange[3],
  Vec3<int32_t> edge1,
  Vec3<int32_t> edge2,
  Vec3<int32_t> v0,
  int poistionClipValue,
  bool haloFlag,
  bool adaptiveHaloFlag,
  bool fineRayflag);

//============================================================================

struct TrisoupSegment {
  Vec3<int32_t> startpos;  // start point of edge segment
  Vec3<int32_t> endpos;    // end point of edge segment

  int index;        // index of segment, to reorder after sorting
  int uniqueIndex;  // index of uniqueSegment
  int vertex;       // distance along segment for intersection (else -1)
};

struct TrisoupSegmentEnc : public TrisoupSegment {
  TrisoupSegmentEnc(
    const Vec3<int32_t>& startpos,
    const Vec3<int32_t>& endpos,
    int index,
    int uniqueIndex,
    int vertex,
    int count,
    int distanceSum,
    int count2,
    int distanceSum2)
    : TrisoupSegment{startpos, endpos, index, uniqueIndex, vertex}
    , count(count)
    , distanceSum(distanceSum)
    , count2(count2)
    , distanceSum2(distanceSum2)
  {}

  int count;        // count of voxels adjacent to this segment
  int distanceSum;  // sum of distances (along segment) of adjacent voxels  
  int count2;        // count of voxels adjacent to this segment
  int distanceSum2;  // sum of distances (along segment) of adjacent voxels  
};

struct TrisoupSegmentNeighbours {
  Vec3<int32_t> startpos;  // start point of edge segment
  Vec3<int32_t> endpos;    // end point of edge segment

  int index;        // index of segment, to reorder after sorting 
  uint16_t neighboursMask;   
};



//----------------------------------------------------------------------------
// comparison for sorting

bool operator<(const TrisoupSegment& s1, const TrisoupSegment& s2);

bool operator<(const TrisoupSegmentNeighbours& s1, const TrisoupSegmentNeighbours& s2);




//============================================================================

}  // namespace pcc
