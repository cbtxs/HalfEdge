
/** 
 * Author : Chunyu Chen 
 * Date   : 2024-4-22
 * Brief  : Base class for traits of halfedge mesh and normal traits.
 */

#include "entity.h"
#include "geometry.h"

namespace HEM
{

/**
 * @brief 半边网格的特性
 */
template<typename N, typename E, typename C, typename H, int D>
class HalfEdgeMeshTraits
{
public:

  using Node = N; 
  using Edge = E; 
  using Cell = C; 
  using HalfEdge = H;

  constexpr static const uint8_t Dim = D;
  static_assert(Dim == 2 || Dim == 3, "Dimension must be 2 or 3.");

  using Point  = typename std::conditional<(D == 2), Point2d, Point3d>::type;
  using Vector = typename std::conditional<(D == 2), Vector2d, Vector3d>::type;
};


/** 
 * @brief 默认 3 维网格的特性
 */
class DefaultNode3d;
class DefaultEdge3d;
class DefaultCell3d;
class DefaultHalfEdge3d;

using DefaultHalfEdgeMesh3dTraits  = HalfEdgeMeshTraits<DefaultNode3d, DefaultEdge3d, DefaultCell3d, DefaultHalfEdge3d, 3>;

class DefaultNode3d : public TNode<DefaultHalfEdgeMesh3dTraits>{};
class DefaultEdge3d : public TEdge<DefaultHalfEdgeMesh3dTraits>{};
class DefaultCell3d : public TCell<DefaultHalfEdgeMesh3dTraits>{};
class DefaultHalfEdge3d : public THalfEdge<DefaultHalfEdgeMesh3dTraits>{};


/** 
 * @brief 默认 2 维网格的特性
 */
class DefaultNode2d;
class DefaultEdge2d;
class DefaultCell2d;
class DefaultHalfEdge2d;

using DefaultHalfEdgeMesh2dTraits  = HalfEdgeMeshTraits<DefaultNode2d, DefaultEdge2d, DefaultCell2d, DefaultHalfEdge2d, 2>;

class DefaultNode2d : public TNode<DefaultHalfEdgeMesh2dTraits>{};
class DefaultEdge2d : public TEdge<DefaultHalfEdgeMesh2dTraits>{};
class DefaultCell2d : public TCell<DefaultHalfEdgeMesh2dTraits>{};
class DefaultHalfEdge2d : public THalfEdge<DefaultHalfEdgeMesh2dTraits>{};

template<int D>
using DefaultHalfEdgeMeshTraits = typename std::conditional<(D == 2), DefaultHalfEdgeMesh2dTraits, DefaultHalfEdgeMesh3dTraits>::type;

}
