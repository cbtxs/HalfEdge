
#include <memory>

#include "Geometry.h"
#include "Vertex.h"
#include "Edge.h"
#include "Cell.h"
#include "HalfEdge.h"
#include "HalfEdgeMeshBase.h"

using Vector = HEM::Vector;
using Vertex = HEM::Vertex;
using Edge = HEM::Edge;
using Cell = HEM::Cell;
using HalfEdge = HEM::HalfEdge;
using HalfEdgeMesh = HEM::HalfEdgeMeshBase<Vertex, Edge, HalfEdge, Cell>;

int main()
{
  auto mesh = std::make_shared<HalfEdgeMesh>();
}
