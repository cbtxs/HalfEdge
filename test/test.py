import ctypes
import numpy as np
import matplotlib.pyplot as plt

from fealpy.mesh.halfedge_mesh import HalfEdgeMesh2d

# Define structures for input and output parameters
class MeshParameter(ctypes.Structure):
    _fields_ = [("a", ctypes.c_double),
                ("b", ctypes.c_double),
                ("c", ctypes.c_double),
                ("d", ctypes.c_double),
                ("nx", ctypes.c_int),
                ("ny", ctypes.c_int)]

class InterfaceParameter(ctypes.Structure):
    _fields_ = [("point", ctypes.POINTER(ctypes.c_double)),
                ("is_fixed_point", ctypes.POINTER(ctypes.c_bool)),
                ("segment", ctypes.POINTER(ctypes.c_int)),
                ("NP", ctypes.c_int),
                ("NS", ctypes.c_int)]

class OutParameter(ctypes.Structure):
    _fields_ = [
                ("point_out1", ctypes.POINTER(ctypes.c_double)),
                ("halfedge_out1", ctypes.POINTER(ctypes.c_int)),
                ("inner_cell1", ctypes.POINTER(ctypes.c_int)),
                ("point_out2", ctypes.POINTER(ctypes.c_double)),
                ("halfedge_out2", ctypes.POINTER(ctypes.c_int)),
                ("idx0", ctypes.POINTER(ctypes.c_int)),
                ("idx1", ctypes.POINTER(ctypes.c_int)),
                ("N", ctypes.POINTER(ctypes.c_int))]


class CutMeshAlgorithm():
    def __init__(self, box, nx, ny, N = -1):
        self.box = box
        self.nx  = nx
        self.ny  = ny
        self.N   = int(1.2*(nx+2)*(ny+2)) if N < 0 else N # 节点的单元个数的上限

        # Load the shared library
        self.lib = ctypes.CDLL("./apps/libcutmesh.so")

        # Define the argument and return types for the C++ functions
        self.lib.get_cut_mesh.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                     ctypes.c_int, ctypes.c_int,
                                     np.ctypeslib.ndpointer(dtype=np.double),
                                     np.ctypeslib.ndpointer(dtype=np.bool_),
                                     np.ctypeslib.ndpointer(dtype=np.intc),
                                     ctypes.c_int, ctypes.c_int,
                                     np.ctypeslib.ndpointer(dtype=np.double),
                                     np.ctypeslib.ndpointer(dtype=np.intc),
                                     np.ctypeslib.ndpointer(dtype=np.intc)]
        self.lib.get_cut_mesh.restype = None

        self.lib.get_cut_mesh2.argtypes = [MeshParameter, 
                                           InterfaceParameter, 
                                           InterfaceParameter, 
                                           OutParameter]
        self.lib.get_cut_mesh2.restype = None


    def get_cut_mesh(self, point, is_fixed_point, segment):
        a, b, c, d = self.box
        nx = self.nx
        ny = self.ny
        
        NP = is_fixed_point.size
        NS = segment.size

        # Prepare the output arrays
        point_out = np.empty(self.N*2, dtype=np.double)
        halfedge_out = np.empty(self.N*24, dtype=np.intc)
        num = np.empty(2, dtype=np.intc)

        # Call the C++ function
        self.lib.get_cut_mesh(a, c, b, d, nx, ny, point, is_fixed_point, segment, 
                NP, NS, point_out, halfedge_out, num)

        # deal data
        node     = point_out[:num[0]].reshape(-1, 2)
        halfedge = halfedge_out[:num[1]].reshape(-1, 6)
        return self._to_halfedge_mesh(node, halfedge)

    def get_cut_mesh2(self, point0, point1, is_fixed_point, segment):
        a, b, c, d = self.box
        nx = self.nx
        ny = self.ny
        
        NP = is_fixed_point.size
        NS = segment.size

        N = self.N

        # Define input parameters
        mesh_param = MeshParameter(a=a, b=c, c=b, d=d, nx=nx, ny=ny)

        #interface_param0 = InterfaceParameter(point0, is_fixed_point, segment, 4, 5)
        interface_param0 = InterfaceParameter(point=point0.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                          is_fixed_point=is_fixed_point.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                                          segment=segment.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                                          NP=NP, NS=NS)
        interface_param1 = InterfaceParameter(point=point1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                          is_fixed_point=is_fixed_point.ctypes.data_as(ctypes.POINTER(ctypes.c_bool)),
                                          segment=segment.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                                          NP=NP, NS=NS)

        point_out1 = np.empty(N*2, dtype=np.double)
        halfedge_out1 = np.empty(N*24, dtype=np.intc)
        inner_cell1 = np.empty(N, dtype=np.intc)

        point_out2 = np.empty(N*2, dtype=np.double)
        halfedge_out2 = np.empty(N*24, dtype=np.intc)

        idx0 = np.empty(N, dtype=np.intc)
        idx1 = np.empty(N, dtype=np.intc)
        num = np.empty(6, dtype=np.intc)

        # Create an instance of OutParameter
        out_param = OutParameter(
            point_out1=point_out1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            halfedge_out1=halfedge_out1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            inner_cell1=inner_cell1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),

            point_out2=point_out2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            halfedge_out2=halfedge_out2.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),

            idx0=idx0.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            idx1=idx1.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            N=num.ctypes.data_as(ctypes.POINTER(ctypes.c_int)))

        # Call the C++ function
        self.lib.get_cut_mesh2(mesh_param, interface_param0, interface_param1, out_param)

        node1 = point_out1[:num[0]].reshape(-1, 2)
        halfedge1 = halfedge_out1[:num[1]].reshape(-1, 6)
        cellflag = inner_cell1[:num[2]]

        node2 = point_out2[:num[3]].reshape(-1, 2)
        halfedge2 = halfedge_out2[:num[4]].reshape(-1, 6)

        mesh1 = self._to_halfedge_mesh(node1, halfedge1)
        mesh2 = self._to_halfedge_mesh(node2, halfedge2)

        idx0 = idx0[:num[5]]
        idx1 = idx1[:num[5]]
        return mesh1, mesh2, idx0, idx1

    def _to_halfedge_mesh(self, node, halfedge):
        _, hidx = np.unique(halfedge[:, 5], return_index=True)
        flag    = halfedge[hidx, 4]!=hidx
        hidx    = np.r_[hidx, halfedge[hidx[flag], 4]]
        hidxmap = np.arange(len(hidx))
        hidxmap[hidx] = np.arange(len(hidx))

        halfedge[:, 2:] = hidxmap[halfedge[:, 2:]] 
        halfedge = halfedge[hidx, :5]

        mesh = HalfEdgeMesh2d(node, halfedge)
        return mesh

# test for get_cut_mesh
point = np.array([0.22131245, 0.21252151, 0.213515125, 0.6566125,
    0.712341251235, 0.65, 0.713515125, 0.23], dtype=np.double)
is_fixed_point = np.array([1, 1, 1, 1], dtype=np.bool_)
segment = np.array([3, 2, 1, 0, 3], dtype=np.intc)

cutalg = CutMeshAlgorithm([0, 1, 0, 1], 10, 10)
mesh = cutalg.get_cut_mesh(point, is_fixed_point, segment)

fig  = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_cell(axes, showindex=True)

# test for get_cut_mesh2
point0 = np.array([0.22131245, 0.21252151, 0.213515125, 0.6566125,
    0.712341251235, 0.65, 0.713515125, 0.23], dtype=np.double)
point1 = point0.copy()
point1[1::2] -= 0.03

mesh1, mesh2, idx0, idx1 = cutalg.get_cut_mesh2(point0, point1, is_fixed_point, segment)

fig  = plt.figure()
axes = fig.gca()
mesh1.add_plot(axes)
mesh1.find_cell(axes, showindex=True)

fig  = plt.figure()
axes = fig.gca()
mesh2.add_plot(axes)
mesh2.find_cell(axes, showindex=True)
plt.show()












